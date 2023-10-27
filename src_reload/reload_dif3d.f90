Subroutine reload_rundif3d(NOUT, num_candidate, candidate_pos, iload, reload_comp, keffs, peakpwr, peakass, &
    &   peakden, peakpos, num_test)
    use m_reload
    use m_utility
    use ZNATDN_io
    implicit none
    integer :: NOUT, num_candidate, iload
    integer :: candidate_pos(num_candidate)
    type(composition_data) :: reload_comp
    real(8) :: keffs(num_candidate), peakpwr(num_candidate), peakden(num_candidate)
    integer :: peakass(num_candidate), peakpos(2,num_candidate), num_test
    ! intrinsic functions
    integer :: chdir  
    ! local
    integer :: isomap(reload_comp%niso), peak_node(num_xynod)
    integer :: ii, jj, IX, IY, IASS
    character(len=8) :: ANAME
    !character(len=100) :: line
    type (ZNATDN_DATA) :: user_ZNATDN
    integer :: niso, ipos, IZ, ireg, izon, last_check
    real(8) :: fdum, max_keff, min_keff, eps
    logical :: check_this, check_next
    integer :: num_accepted
    character(len=16) :: objective

    100 format('[reload_rundif3d]...',A,1X,I6)

    if(.not. NDXSRF%defined) then     ! first entry
        ! load base NDXSRF and ZNATDN
        call NDXSRF_IMPORT(NOUT, NDXSRF, 'NDXSRF')
        call ZNATDN_IMPORT(NOUT, ZNATDN, 'ZNATDN')  ! original ZNATDN
        ! modify NDXSRF, ZNATDN, and LABELS (slightly) by removing subzones and using a single nuclide set
        if(NDXSRF%NON > NDXSRF%NNS .or. NDXSRF%NSN > 0 .or. NDXSRF%NSZ > 0) &
            call ARC_Simplify_NDXSRF(NOUT, NDXSRF, ZNATDN, LABELS)
        call NDXSRF_EXPORT(NOUT, NDXSRF, 'user_NDXSRF') ! in case it is modified
        call LABELS_EXPORT(NOUT, LABELS, 'user_LABELS') ! in case it is modified
    end if
    call system_command(NOUT, 'cp ISOTXS '//dif3d_dir//'/ISOTXS')
    call system_command(NOUT, 'cp GEODST '//dif3d_dir//'/GEODST')
    call system_command(NOUT, 'cp user_NDXSRF '//dif3d_dir//'/NDXSRF')
    call system_command(NOUT, 'cp user_LABELS '//dif3d_dir//'/LABELS')

    ! move in DIF3D temporary working directory
    ii = chdir(dif3d_dir)
    if(ii /= 0) then
        write(NOUT,100) 'Error: failed to change working directory to '//dif3d_dir
        call abort 
    end if
    
    ! locate nuclide in NDXSRF nuclide set 
    niso = reload_comp%niso
    do ii = 1, niso
        ANAME = reload_comp%isonam(ii)
        do jj = 1, NDXSRF%NON
            if(ANAME == NDXSRF%HNNAME(jj)) exit 
        end do
        if(jj > NDXSRF%NON) then
            write(NOUT,100) 'Error: isotope '//ANAME//' cannot be found in NDXSRF'
            call abort
        end if
        isomap(ii) = jj
    end do
    
    ! loop over candidate reloading positions
    if(consider_keff) then
        max_keff = eoc_keff + max_deltak
        min_keff = eoc_keff + min_deltak
    end if
    keffs(:) = -1.0
    num_accepted = 0
    last_check = candidate_pos(1)
    do ipos = 1, num_candidate
        !write(0,*) 'working on candidate',ipos,' assembly ',candidate_pos(ipos)
        check_this = .true.
        ! if this candidate position gives practically same estimation as last one, skip expensive DIF3D calculation for it
        if((.not. dif3d_check_all) .and. ipos > 1 .and. ((objective_func /= 'reactivity' .and. objective_func /= 'reactivity_slow') .or. use_persent)) then
            if(objective_func == 'assembly_power') then
                fdum = pekpwr(candidate_pos(ipos),iload) / pekpwr(last_check,iload) - 1.0D0
                eps = eps_power
                objective = 'assembly power'
            else if(objective_func == 'peak_power') then
                fdum = pekden(candidate_pos(ipos),iload) / pekden(last_check,iload) - 1.0D0
                eps = eps_power
                objective = 'peak power'
            else if(objective_func == 'burnup') then
                fdum = assbrn(candidate_pos(ipos)) / assbrn(last_check) - 1.0D0
                eps = eps_power
                objective = 'burnup'
            else if(use_persent) then
                fdum = deltak(candidate_pos(ipos),iload) - deltak(last_check,iload)
                eps = eps_keff
                objective = 'delta_k'
            end if
            if(abs(fdum) <= eps) then  ! this position will probably give similar peak powers
                peakpwr(ipos) = peakpwr(ipos-1)
                peakass(ipos) = peakass(ipos-1)
                peakden(ipos) = peakden(ipos-1)
                peakpos(:,ipos) = peakpos(:,ipos-1)
                keffs(ipos) = keffs(ipos-1)
                write(NOUT,'(" skipped DIF3D calculation for assembly ",I4," because it gives practically same estimated ",A, &
                    &" as the last candidate, assembly ",I4)') candidate_pos(ipos), trim(objective), last_check
                check_this = .false.
            else
                last_check = abs(candidate_pos(ipos))
            end if
        end if

        if(check_this) then
            !write(0,*) '    run dif3d for this candidate'
            IASS = candidate_pos(ipos)
            call reload_ass2xy(IASS, num_xnode, IX, IY)
            ! modify ZNATDN for reloaded regions, assuming that different assemblies will not have same fuel zone (almost always true for REBUS model)
            ! otherwise, need to add a new zone for reloaded assembly and modify GEODST, NDXSRF, and LABELS accordingly (not implemented yet)
            call ZNATDN_Copy(ZNATDN, user_ZNATDN)
            do IZ = bot_active_node, top_active_node
                ireg = GEODST%MR(IX,IY,IZ)
                if(IZ /= 1 .and. ireg == GEODST%MR(IX,IY,IZ-1)) cycle   ! same region as last node
                izon = GEODST%NZNR(ireg)    ! one-to-one correspondence between fuel region and zone in REBUS
                user_ZNATDN%ADEN(:,izon) = 0.0E0
                do ii = 1, niso
                    jj = isomap(ii)
                    user_ZNATDN%ADEN(jj,izon) = reload_comp%density(ii)
                end do
            end do
            ! run DIF3D
            call ZNATDN_EXPORT(NOUT, user_ZNATDN, 'ZNATDN')
            call system_command(NOUT, '../dif3d.x < dif3d.inp > dif3d.out')
            ! get keff and peak assembly power
            call read_boc_keff('dif3d.out', NOUT, keffs(ipos))
            call PWDINT_IMPORT_BINARY(NOUT, boc_PWDINT, 'PWDINT')
            call get_peak_power(boc_PWDINT,NOUT,bocpwr,peakass(ipos),peakpwr(ipos),bocden,peak_node,peakpos(:,ipos),peakden(ipos))
        end if

        check_next = .false.
        ! check whether all constraints are satisfied
        fdum = peakden(ipos) / peak_power_limit - 1.0D0
        if(fdum > eps_power) then
            if(fixed_peak_power .or. objective_func /= 'peak_power') check_next = .true.
        end if
        if(consider_keff) then
            if(keffs(ipos) - max_keff > eps_keff) check_next = .true.
            if(mindk > 0.0 .and. keffs(ipos) - min_keff < -eps_keff) check_next = .true.
        end if
        if(.not. check_next) then
            num_accepted = num_accepted + 1
        else 
            candidate_pos(ipos) = -candidate_pos(ipos)  ! mark candidate not meet all constraints
        end if
        if(ipos == num_candidate) check_next = .false.
        if(dif3d_check_all .and. ipos < num_candidate) cycle

        ! extend checking for more candidates under constraint:
        ! giving 0~10% higher estimated peak powers than minimum
        if(ipos < num_candidate .and. (objective_func=='assembly_power' .or. objective_func=='peak_power')) then
            if(objective_func == 'assembly_power') then
                fdum = pekpwr(candidate_pos(ipos+1),iload) / pekpwr(abs(candidate_pos(1)),iload) - 1.0D0 
            else
                fdum = pekden(candidate_pos(ipos+1),iload) / pekden(abs(candidate_pos(1)),iload) - 1.0D0 
            end if
            if(fdum <= power_check_range .or. num_accepted == 0) then
                check_next = .true.
            else ! (fdum > power_check_range .and. num_accepted > 0)
                check_next = .false.
            end if
        end if
        ! giving 0~10% lower burnup than maximum
        if(objective_func == 'burnup' .and. ipos < num_candidate) then
            fdum = assbrn(candidate_pos(ipos+1)) / assbrn(abs(candidate_pos(1))) - 1.0D0
            if(fdum >= -burnup_check_range .or. num_accepted == 0) then
                check_next = .true.
            else
                check_next = .false.
            end if
        end if

        ! extending check range to candidates with 0~min_deltak lower reactivity 
        if(objective_func == 'reactivity' .and. ipos < num_candidate) then
            if(use_persent) then
                fdum = deltak(abs(candidate_pos(1)),iload) - deltak(candidate_pos(ipos+1),iload)
                if(fdum - min_deltak < -eps_keff .or. num_accepted == 0) then
                    check_next = .true.
                else
                    check_next = .false.
                end if
            else        ! no idea how to check without estimated keff
                if(num_accepted < 3) check_next = .true.
            end if
        end if

        if(objective_func == 'reactivity_slow' .and. ipos < num_candidate) then
            if(use_persent) then
                fdum = deltak(abs(candidate_pos(1)),iload) - deltak(candidate_pos(ipos+1),iload) - min_deltak             ! 
                fdum = fdum - conv_peak_rho * wgt_pekfac * (pekfac(abs(candidate_pos(1)),iload) - pekfac(candidate_pos(ipos+1),iload)) ! conversion 1 pcm reactivity = 1 percent peaking factor
                if(fdum < -eps_keff .or. num_accepted == 0) then
                    check_next = .true.
                else
                    check_next = .false.
                end if
            else        ! no idea how to check without estimated keff
                if(num_accepted < 3) check_next = .true.
            end if
        end if
        
        if(check_next) cycle
        ! sucessful return with a verified refueling position
        num_test = ipos
        exit

        ! temporary use: save dif3d output
        !write(line,'(A,I3.3)') 'mv dif3d.out ../dif3d.out_',IASS
        !call system_command(NOUT, trim(line))
    end do

    ! now move back to rebus directory
    ii = chdir('..')
    if(ii /= 0) then
        write(NOUT,100) 'Error: failed to move back to REBUS working directory'
        call abort 
    end if

End Subroutine


! ----------------------------------------------------------------------
Subroutine reload_dif3d_input(NOUT, filename, init)
    ! make dif3d input with existing interface datasets
    use m_utility
    use NE_Kind, only: NE_Kind_GetFreeLogicalUnit
    implicit none
    integer :: NOUT
    integer :: init   ! flag for initial guess, init = 0 means no initial flux guess
    !   init =  1 / 2 / 3  initial flux guess exists as NHFLUX / NAFLUX / both    
    !   existing interface files including ISOTXS, NDXSRF, ZNATDN, GEODST, LABELS
    character(len=*) :: filename
    ! local
    integer :: UNIT, NIN, IDUM, freeform
    character(len=100) :: line
    character(len=8) :: ANAME, datasets(5), files(2)
    integer :: ncards(99), MCARD, ii, jj
    logical :: yes
    data datasets /'ISOTXS','GEODST','LABELS','NDXSRF','ZNATDN'/
    data files /'ADIF3D', 'AHMG4C'/
    
    100 format('[reload_dif3d_input]...',A,1X,I6)
    
    UNIT = NE_Kind_GetFreeLogicalUnit()
    ! make dif3d input with existing interface datasets
    open(UNIT,file=filename,status='replace',form='formatted')
    write(UNIT,'("BLOCK=OLD")')
    do ii = 1, 5
        write(UNIT,'("DATASET=",A)') datasets(ii)
    end do
    if(init == 1 .or. init == 3) write(UNIT,'("DATASET=NHFLUX")')
    if(init == 2 .or. init == 3) write(UNIT,'("DATASET=NAFLUX")')
    write(UNIT,'("BLOCK=STP021,3")')

    NIN = NE_Kind_GetFreeLogicalUnit()
    do jj = 1, 2
        inquire(file=files(jj), exist=yes)
        if(yes) then
            open(NIN, file=files(jj), status='old', form='formatted')
            ! copy input deck
            call scan_arc_file_header(NIN,MCARD,ncards,freeform,ANAME)
            if(freeform==1) then
                write(UNIT,'("UNFORM=",A)') trim(ANAME)
            else
                write(NOUT,100) 'Error: '//ANAME//' of DIF3D input has to be in free format to be used by PERSENT'
                call abort
            end if
            MCARD = sum(ncards(1:MCARD))
            do ii = 1, MCARD
                read(NIN,'(A)') line
                read(line(1:2),*) IDUM
                if(jj==1 .and. IDUM == 4) then    ! force to export PWDINT for obtaining peak node power
                    write(UNIT,'("04         0     0     0    00   000    10   000     0     0  0001     0")')
                else
                    write(UNIT,'(A)') trim(line)
                end if
            end do
            if(jj==1 .and. ncards(4)==0) write(UNIT,'("04         0     0     0    00   000    10   000     0     0  0001     0")')
            close(NIN)
        else if(jj==2) then  ! A.HMG4C is optional
            write(UNIT,'("UNFORM=A.HMG4C",/,"01    DISABLE HMG4C EDITS",/,&
                "02  20000000   1     0     0     0     1",/)')
        else
            write(NOUT,100) 'Error: cannot find DIF3D interface file ADIF3D'
            call abort
        end if
    end do
    close(UNIT)

End Subroutine

! ----------------------------------------------------------------------
