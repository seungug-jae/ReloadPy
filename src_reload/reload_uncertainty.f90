subroutine reload_uncertainty(NOUT, reload_pos)
    use m_reload
    use m_utility
    use NE_Kind, only: NE_Kind_GetFreeLogicalUnit
    implicit none
    ! character(len=12), parameter :: persent_dir   = 'work_persent'
    integer :: NOUT, reload_pos
    ! local
    integer :: ii, jj, kk, IASS, IX, IY,iload
    integer :: NIN, MCARD, ncards(99), freeform, HMN
    integer :: chdir    ! intrinsic function
    character(len=8) :: ANAME
    logical :: forward_avaiable = .false.
    logical :: adjoint_avaiable = .false.
    logical :: adjoint_initial = .false.
    real(8) :: fdum, fdum2, discharged, refueled
    integer :: num_cand_, cand_pos_(1) ! dummy data for executing files
    real(8), allocatable :: deltak_in(:,:)    ! (num_xynod, num_reload_comp+1) change in keff due to composition replacement of specific fuel assembly
    real(8), allocatable :: asswth_in(:,:)    ! (num_xynod, num_reload_comp) total assembly worth of reloaded fuel assembly
    real(8), allocatable :: asspwr_in(:,:)   ! (num_xynod, num_reload_comp) powers of reloaded assembly at different positions
    real(8), allocatable :: pekpwr_in(:,:)   ! (num_xynod, num_reload_comp) peak assembly powers with reloaded assembly at different positions
    integer, allocatable :: pekass_in(:,:)   ! (num_xynod, num_reload_comp) position (GEODST ordering) of assembly having peak power with reloaded assembly
    real(8), allocatable :: pekden_in(:,:)   ! (num_xynod, num_reload_comp) peak (node-averaged) power density corresponding to different reloading positions
    integer, allocatable :: pekpos_in(:,:,:) ! (2,num_xynod, num_reload_comp) node position (IASS,IZ) showing peak power density corresponding to different reload positions
    real(8) :: peaks(num_potential_site)   ! used as scratch working array at first and then store peak assembly powers obtained with DIF3D calculations
    real(8) :: maxup(num_potential_site)   ! maximum incurease in assembly power of reloaded core from reference condition
    real(8) :: rmsdif(num_potential_site)  ! root mean square change in assembly power distribution
    integer :: maxpos
    real(8) :: avgpwr, digit_shift

    100 format('[reload_uncertainty]...',A,1X,I6)

    num_cand_ = 1
    cand_pos_(1) = reload_pos

    !-----------------------------------------------------------------------
    if (consider_worth) then
        IF (allocated(deltak_in)) DEALLOCATE(deltak_in)
        IF (allocated(asswth_in)) DEALLOCATE(asswth_in)
        allocate(deltak_in(num_xynod, num_test_comp))
        allocate(asswth_in(num_xynod, num_test_comp))
        deltak_in(:,:) = 0.0d0
        asswth_in(:,:) = 0.0d0
    end if

    if(use_persent) then
        ! check VARIANT option
        inquire(file='user.NHFLUX', exist=forward_avaiable)
        inquire(file='user.NAFLUX', exist=adjoint_avaiable)
        if(.not. forward_avaiable) then
            NIN = NE_Kind_GetFreeLogicalUnit()
            open(NIN,file='ADIF3D',status='old',form='formatted')
            call scan_arc_file_header(NIN,MCARD,ncards,freeform,ANAME)
            if(freeform /= 1) then
                write(NOUT,100) 'Error: when PERSENT is invoked, A.DIF3D input has to be in free format'
                call abort 
            end if
            if(MCARD >= 12 .and. ncards(12) > 0) then 
                ii = sum(ncards(1:11))
                call SkipLines(NIN, ii)
                read(NIN,*) ii, ii, HMN
                HMN = abs(HMN)
                if(HMN > 9999) HMN = modulo(HMN,10000)
                if(HMN > 99) then
                    ii = HMN / 100              ! M = N is guaranteed 
                else
                    ii = HMN / 10
                end if
                if(ii == 1) then
                    forward_avaiable = .true.     ! see below, if NHFLUX at EOC is saved, we can fix this 
                    call system_command(NOUT,'ln -s  NHFLX0  user.NHFLUX')
                end if
            end if
            close(NIN)
        end if

        ! set up persent run environment
        if(adjoint_avaiable) then
            call system_command(NOUT, 'rm '//persent_dir//'/*')    ! not first run
        else
            call system_command(NOUT, 'mkdir '//persent_dir)
            inquire(file='init.NAFLUX', exist=adjoint_initial)
            if(adjoint_initial) call system_command(NOUT, 'cp init.NAFLUX '//persent_dir//'/NAFLUX')   ! used as initial guess for adjoint calculation
        end if
        call system_command(NOUT, 'cp ISOTXS '//persent_dir//'/ISOTXS')         ! for DIF3D null run
        call system_command(NOUT, 'cp ISOTXS '//persent_dir//'/user.ISOTXS')    ! for PERSENT
        call system_command(NOUT, 'cp GEODST '//persent_dir//'/GEODST')
        call system_command(NOUT, 'cp NDXSRF '//persent_dir//'/NDXSRF')
        call system_command(NOUT, 'cp LABELS '//persent_dir//'/LABELS')
        call system_command(NOUT, 'cp ZNATDN '//persent_dir//'/ZNATDN')

        ! try to utilize REBUS output, but it seems REBUS has problem with saving NHFLUX at EOC. Only NHFLX0 is updated at EOC
        if(forward_avaiable) then
            if(adjoint_avaiable) then
                call reload_persent_input_in('persent_tcomp.inp',2, num_cand_, cand_pos_)
            else
                call reload_persent_input_in('persent_tcomp.inp',1, num_cand_, cand_pos_)
            end if
        else
            call reload_persent_input_in('persent_tcomp.inp',0, num_cand_, cand_pos_)   ! let PERSENT generate NHFLUX and NAFLUX
        end if
        ! ! make dif3d input using interface files
        if(adjoint_initial) then
            call reload_dif3d_input(NOUT,'reload.dif3d',2)
        else
            call reload_dif3d_input(NOUT,'reload.dif3d',0)
        end if
        call system_command(NOUT, 'cp reload.dif3d '//persent_dir//'/dif3d.inp')
        call system_command(NOUT, 'cp persent_tcomp.inp '//persent_dir//'/persent_tcomp.inp')
        
        ! move in temporary working directory
        ii = chdir(persent_dir)
        if(ii /= 0) then
            write(NOUT,100) 'Error: failed to change working directory to '//persent_dir
            call abort 
        end if

        ! link existing flux data (avoid copy of large files)
        if(forward_avaiable) call system_command(NOUT, 'ln -s ../user.NHFLUX .')
        if(adjoint_avaiable) call system_command(NOUT, 'ln -s ../user.NAFLUX .')

        ! run persent
        call CPU_TIME(fdum)
        call system_command(NOUT, '../persent.x  persent_tcomp.inp > persent_tcomp.out')
        call CPU_TIME(fdum2)
        persent_time(2) = persent_time(2) + fdum2 - fdum
        ! grab reactivity perturbation results for each fuel assembly from persent output
        call CPU_TIME(fdum)
        call reload_persent_output_in(NOUT, 'persent_tcomp.out', eoc_keff)
        call CPU_TIME(fdum2)
        persent_time(3) = persent_time(3) + fdum2 - fdum

        ! ! save NHFLUX and NAFLUX in case additional PERSENT run is needed
        if(.not. forward_avaiable) call system_command(NOUT, 'mv base.NHFLUX ../user.NHFLUX')
        if(.not. adjoint_avaiable) call system_command(NOUT, 'mv base.NAFLUX ../user.NAFLUX')
        
        ! compute assembly worth 
        if(consider_worth) then
            jj = num_reload_comp + 1
            do ii = 1, num_test_comp
                ! do kk = 1, num_cand_
                IASS = reload_pos
                asswth_in(IASS,ii) = deltak_in(IASS,ii) - deltak(IASS,jj)
                ! end do
            end do
        end if

        ! move back to REBUS working directory
        ii = chdir('..')
        if(ii /= 0) then
            write(NOUT,100) 'Error: failed to move back to REBUS working directory'
            call abort 
        end if
    end if    
    ! power calculation ------------------------------------------------------------------------------------------------
    if(allocated(asspwr_in)) deallocate(asspwr_in)
    if(allocated(pekpwr_in)) deallocate(pekpwr_in)
    if(allocated(pekass_in)) deallocate(pekass_in)
    if(allocated(pekden_in)) deallocate(pekden_in)
    if(allocated(pekpos_in)) deallocate(pekpos_in)
    
    ! ! if (measure_uncertainty) then
    allocate(asspwr_in(num_xynod,num_test_comp), &
        pekpwr_in(num_xynod,num_test_comp), pekass_in(num_xynod,num_test_comp), &
        pekden_in(num_xynod,num_test_comp), pekpos_in(2,num_xynod,num_test_comp), stat=ii)
    call reload_power_in(NOUT, reload_pos)

    ! print table header
    write(NOUT,'(/,20X,"*** Estimating reactivity change due to reloaded composition uncertainty (upto 4-sigma)   ***")')
    
    if(num_hex_sector > 0) then
        write(NOUT, '(A, I3, 5x, A, I6,2I5)') '   num_test_comp : ', num_test_comp, '   Reloaded position : ', IASS, hexagon_rp(:,IASS)
    else
        call reload_ass2xy(IASS, num_xnode, IX, IY)
        write(NOUT, '(A, I3, 5x, A, I6,2I5)') '   num_test_comp : ', num_test_comp, '   Reloaded position : ', IASS, IX, IY
    end if

    ! if(consider_worth) then
    !     write(NOUT,'(" Note: assembly worth is estimated as k-effective change between two composition replacement with ",&
    !     &   A8," and ",A8,/)') TFUEL(iload)%name, TFUEL(num_reload_comp+1)%name
    ! else
    !     write(NOUT,'(" Note: assembly worth is not calculated (only net change in keff is determined)",/)')
    ! end if
    if(num_hex_sector > 0) then
        write(NOUT,'("          Assembly worth at this position       Added worth    RMS change       Max. increase    Assembly seeing  | Peak Assembly  Corresponding peak",/, &
        &            "MATERIAL  Burned assembly   Fresh assembly      by fueling     /Avg power (%)   /Avg power (%)   maximum increase | Position       assembly power (W)",/)')
    else
        write(NOUT,'("          Assembly worth at this position       Added worth    RMS change       Max. increase    Assembly seeing  | Peak Assembly  Corresponding peak",/, &
        &            "MATERIAL  Burned assembly   Fresh assembly      by fueling     /Avg power (%)   /Avg power (%)   maximum increase | Position       assembly power (W)",/)')
    end if
    ! write(NOUT,'("candidate reload position   RMS change       Max. increase    Assembly seeing    Corresponding peak")')
    ! write(NOUT,'("   IASS    IRING   IPOS     /Avg power (%)   /Avg power (%)   maximum increase   assembly power (W)")')
    avgpwr = PWDINT%POWER / num_potential_site  ! approximate average fuel assembly power
    digit_shift = 10**(-floor(log10(eps_power_diff)))      ! reasonable assumption: eps_power_diff < 1.0
    IASS = reload_pos
    call reload_ass2xy(IASS, num_xnode, IX, IY)
    DO iload = 1, num_test_comp
        call reload_power_change_in(iload, IASS, maxup(1), maxpos, rmsdif(1))
        if(consider_worth) then
            discharged = asswth_in(IASS,iload) - deltak_in(IASS,iload)      ! asswth is worth of fresh fuel assembly if reloaded, deltak is added worth
            refueled = asswth_in(IASS,iload)
        else
            discharged = 0.D0
            refueled = 0.D0            
        end if
    ! ! since many candidates may have similar rmsdif, sorting first by increasing rmsdif and then by increasing peak assembly power
        fdum = pekpwr_in(IASS,iload) / avgpwr      ! should be well below 10.0
        peaks(1) = floor(rmsdif(1)*digit_shift+0.5d0)*10 + fdum
    ! edit
        if(num_hex_sector > 0) then
            write(NOUT,'(A6,6X,F13.5,5X,F13.5,4X,F13.5,5X,2(F10.2,7X),I10,9X,I8, 4X,ES13.5)') &
            TFUEL(iload)%name, discharged,refueled,deltak_in(IASS,iload), rmsdif(1)*100, &
            maxup(1)*100, maxpos, pekass_in(IASS,iload), pekpwr_in(IASS,iload)
        else
            write(NOUT,'(A6,6X,F13.5,5X,F13.5,4X,F13.5,5X,2(F10.2,7X),I10,9X,I8, 4X,ES13.5)') &
            TFUEL(iload)%name, discharged,refueled,deltak_in(IASS,iload), rmsdif(1)*100, &
            maxup(1)*100, maxpos, pekass_in(IASS,iload), pekpwr_in(IASS,iload)
        end if
        ! write(NOUT,'(I7,2X,2I7,5X,2(F10.2,7X),I10,9X,ES13.5)') IASS, hexagon_rp(:,IASS), rmsdif(1)*100, &
        !     maxup(1)*100, maxpos, pekpwr_in(IASS,iload)!, peaks(1)
    END DO

    Contains

    Subroutine reload_persent_output_in(NOUT, filename, user_keff)
        ! read persent output for region wise reactivity worth due to fuel reloading
        use m_reload
        use m_utility
        use NE_Kind, only: NE_Kind_GetFreeLogicalUnit
        implicit none
        integer :: NOUT
        character(len=*) :: filename
        real(8) :: user_keff
        ! local
        integer :: UNIT, IOS, ii, jj, IASS, iload
        real(8) :: drho, keff
        character(len=8) :: ANAME
        integer, parameter :: line_len = 300
        character(len=line_len) :: line
        character(len=1) :: splitter
        
        100 format('[reload_persent_output_in/reload_uncertainty]...',A,1X,I6)
    
        ! read PERSENT output for delta_rho and accumulate reactivity worth by replacing each fuel assembly
        UNIT = NE_Kind_GetFreeLogicalUnit()
        open(UNIT,file=filename, status='old',form='formatted')
    
        keff = user_keff    ! if keff < 0.0,  will get it from PERSENT output
        do jj = 1, num_test_comp
            call SkipTo('[PERSENT]...Problem', .true., UNIT, IOS)
            call SkipTo('[PERSENT]...Performing the DIF3D-VARIANT numerator/denominator operations',.true., UNIT,IOS)
            read(UNIT,'(A)') line
            !ii = index(line,'is')
            !read(line(ii+2:line_len),*) drho
            ! seems fixed format is used in PERSENT output
            ii = index(line,'RF_FA_')
            read(line(ii+6:line_len),'(I5,A1,I2)') IASS, splitter, iload
            read(line(58:line_len),*) drho
            if(keff < 0.0) read(line(120:line_len),*) keff
            deltak_in(IASS,iload) = keff * keff * drho / (1.0d0 - keff*drho)
        end do
        close(UNIT)
    
    End Subroutine
    
    
    ! ----------------------------------------------------------------------
    Subroutine reload_persent_input_in(filename, iflux, num_pool, pool)
        ! make persent input for zone perturbation problem
        use m_reload
        use NE_Kind, only: NE_Kind_GetFreeLogicalUnit
        implicit none
        integer :: iflux      ! iflux = 0 / 1 / 2  prepare input with no flux / forward flux / both forward and adjoint flux provided
        character(len=*) :: filename
        integer :: num_pool, pool(num_pool)
        ! local
        integer :: UNIT, NIN, IDUM, IOS, NCMP
        integer :: ii, jj, kk, ij, iload, IASS, IX, IY, IZ, NZ, IREG, num_pert_region
        character(len=8) :: ANAME
        character(len=8), allocatable :: pert_regions(:)   ! fuel regions in an assembly
        
        UNIT = NE_Kind_GetFreeLogicalUnit()
        open(UNIT,file=filename,status='replace',form='formatted')
        
        if(iflux /= 0) write(UNIT,'("FORCE_FULL_FLUX   NO")')
        !write(UNIT,'("FORCE_FULL_FLUX   NO")')
        write(UNIT,'("DIF3D_EXECUTABLE  ../dif3d.x")')
        write(UNIT,'("DIF3D_INPUT       dif3d.inp")') 
        write(UNIT,'("ISOTXS_INPUT      user.ISOTXS")')
        if(iflux > 0) write(UNIT,'("FORWARD_FILE      user.NHFLUX")')
        if(iflux > 1) write(UNIT,'("ADJOINT_FILE      user.NAFLUX")')
        write(UNIT,*)

        NCMP = num_test_comp
        do iload = 1, NCMP
            ANAME = TFUEL(iload)%name
            jj = (TFUEL(iload)%niso + 9) / 10
            kk = 0
            do ii = 1, jj-1
                write(UNIT,'("NEW_ZONE  ",A,1X,10(1X,A8,ES13.5))') ANAME, ((TFUEL(iload)%isonam(kk+ij), TFUEL(iload)%density(kk+ij)),ij=1,10)
                kk = kk + 10
            end do
            if(kk < TFUEL(iload)%niso) &
            write(UNIT,'("NEW_ZONE  ",A,1X,10(1X,A8,ES13.5))') &
            &   ANAME, ((TFUEL(iload)%isonam(ij), TFUEL(iload)%density(ij)),ij=kk+1,TFUEL(iload)%niso)
            write(UNIT,*)
        end do
        ! end if
        NZ = top_active_node - bot_active_node + 1
        allocate(pert_regions(NZ))
    
        do jj = 1, num_pool
            IASS = pool(jj)
            call reload_ass2xy(IASS, num_xnode, IX, IY)
            num_pert_region = 0
            do IZ = bot_active_node, top_active_node
                IREG = GEODST%MR(IX,IY,IZ)
                if(IZ > bot_active_node .and. IREG == GEODST%MR(IX,IY,IZ-1)) cycle
                ANAME = LABELS%REGNAM(IREG)
                do ii = 1, num_fuel_region
                    if(ANAME == fuel_regions(ii)) exit
                end do
                if(ii <= num_fuel_region) then
                    num_pert_region = num_pert_region + 1
                    pert_regions(num_pert_region) = ANAME
                end if
            end do        
            do iload = 1, NCMP
                ANAME = TFUEL(iload)%name
                do ii = 1, num_pert_region
                    write(UNIT,'("ADJUST_ZONE   RF_FA_",I0.5,"_",I0.2,2X,"FIRST_ORDER_PT",2X,A8,1X,A8," 1.0")') &
                    &   IASS, iload, pert_regions(ii), ANAME
                end do
                write(UNIT,*)
            end do
        end do
        ! num_pert_cases = num_pool * NCMP   ! for the time being, each assembly has at most 2 perturbations, one for fresh fuel and one for empty
    
        close(UNIT)
    
    End Subroutine
    
    Subroutine reload_power_in(NOUT,reload_pos)
        use m_reload
        use NE_Kind, only: NE_Kind_GetFreeLogicalUnit
        implicit none
        integer :: NOUT, reload_pos
        ! local
        integer :: NINTI, NINTJ, NINTK, UNIT
        integer :: IASS, IZ, IX, IY, ii, jj, iload
        real(8) :: AREA, normfac, FDUM
        real(8), parameter :: sqrt3 = 1.732050807568877D0
        logical :: refpwr_exist, write_refpwr
        integer :: peak_ass, peak_pos(2), peak_ass2, peak_pos2(2)
        integer :: local_peak_nod, new_peak_ass, new_peak_pos(2)
        real(8) :: peak_pwr, peak_pwr2, peak_den, peak_den2 ! peak assembly power (W), peak power density (W/cc)
        real(8) :: local_peak_den, new_peak_pwr, new_peak_den
    
        100 format('[reload_power_in/reload_uncertainty]...',A,1X,I6)
        NINTI = GEODST%NINTI
        NINTJ = GEODST%NINTJ
        NINTK = GEODST%NINTK
    
        ! get EOC peak powers
        call get_peak_power(PWDINT, NOUT, eocpwr, peak_ass, peak_pwr, eocden, peknod, peak_pos, peak_den)
        ! edit EOC powers and find the second largest assembly/node powers 
        peak_pwr2 = 0.0D0
        peak_den2 = 0.0D0
        IASS = 0
        ! PRINT *, num_potential_site, peak_ass, peak_pos, eocpwr(peak_ass)
        do ii = 1, num_potential_site
            IASS = fuel_pos(ii)
            if(eocpwr(IASS) >= peak_pwr .and. IASS /= peak_ass) then
                peak_pwr2 = eocpwr(IASS)             ! second largest assembly power
                peak_ass2 = IASS
            end if
            if(eocden(IASS) >= peak_den .and. IASS /= peak_pos(1)) then 
                peak_den2 = eocden(IASS)             ! second largest power density in a different assembly
                peak_pos2(1) = IASS
                peak_pos2(2) = peknod(IASS)
            end if
        end do
    
        normfac = 1.0D0
        if(normalize_power) then
            write(NOUT,'("Renormalize assembly powers after refueling? Yes",/)')
        else
            write(NOUT,'("Renormalize assembly powers after refueling? No",/)')
        end if
        do iload = 1, num_test_comp
            if(flux_data == 'RTFLUX') then
                ! do IY = 1, NINTJ 
                !     do IX = 1, NINTI
                !IASS = (IY-1)*NINTI + IX
                IASS = reload_pos
                call reload_ass2xy(IASS, num_xnode, IX, IY)
                asspwr_in(IASS,ILOAD) = 0.0D0
                pekpwr_in(IASS,ILOAD) = 0.0D0
                pekden_in(IASS,ILOAD) = 0.0D0
                local_peak_den = 0.0D0
                ! if(assembly_map(IX, IY) /= reload_type) cycle
                ! compute power for reloaded assembly under same flux
                asspwr_in(IASS,ILOAD) = eocpwr(IASS)
                do IZ = bot_active_node, top_active_node
                    FDUM = sum(kerma%nheat(:,num_reload_comp+iload) * RTFLUX%FREG(IX,IY,IZ,:))
                    asspwr_in(IASS,ILOAD) = asspwr_in(IASS,ILOAD) + (FDUM - PWDINT%PWR(IX,IY,IZ)) * node_vol(IZ)
                    if(FDUM > local_peak_den) then
                        local_peak_den = FDUM
                        local_peak_nod = IZ
                    end if
                end do
                if(normalize_power) then
                    normfac = PWDINT%POWER / (PWDINT%POWER - eocpwr(IASS) + asspwr_in(IASS,ILOAD))   ! renormalization factor
                    !write(NOUT,'("reloaded at",I6,"  renomalization factor =",ES13.5)') IASS, normfac
                end if
                ! IF (iass==1262) print *, normfac, asspwr_in(IASS,iload)
                asspwr_in(IASS,ILOAD) = asspwr_in(IASS,ILOAD) * normfac
                local_peak_den = local_peak_den * normfac
                if(IASS == peak_ass) then
                    new_peak_pwr = peak_pwr2 * normfac
                    new_peak_ass = peak_ass2
                else
                    new_peak_pwr = peak_pwr * normfac
                    new_peak_ass = peak_ass
                end if
                if(IASS == peak_pos(1)) then
                    new_peak_den = peak_den2 * normfac 
                    new_peak_pos = peak_pos2
                else
                    new_peak_den = peak_den * normfac
                    new_peak_pos = peak_pos
                end if
                if(asspwr_in(IASS,ILOAD) > new_peak_pwr) then
                    new_peak_pwr = asspwr_in(IASS, iload)
                    pekass_in(IASS,ILOAD) = IASS
                else
                    pekass_in(IASS,ILOAD) = new_peak_ass
                end if
                if(local_peak_den > new_peak_den) then
                    new_peak_den = local_peak_den
                    pekpos_in(:,IASS,iload) = (/IASS, local_peak_nod/)
                else
                    pekpos_in(:,IASS,iload) = new_peak_pos
                end if
                pekpwr_in(IASS,ILOAD) = new_peak_pwr       ! peak assembly-integrated power
                pekden_in(IASS,ILOAD) = new_peak_den       ! peak node-averaged power density
                !     end do
                ! end do
            else if(flux_data == 'NHFLUX') then
                write(NOUT,100) 'Error: usage of NHFLUX not implemented yet'
                call abort
            end if
        end do
    
    End Subroutine
    Subroutine reload_power_change_in(iload, reload_pos, maxup, maxpos, rmsdif)
        use m_reload
        implicit none
        integer :: iload, reload_pos, maxpos
        real(8) :: maxup, rmsdif  ! max increase, root mean square difference
        ! local
        real(8) :: pwrdif, scale, inverse_avgpwr
        integer :: ii, IASS
    
        scale = (PWDINT%POWER - asspwr_in(reload_pos,iload)) / (PWDINT%POWER - eocpwr(reload_pos))
        inverse_avgpwr = num_potential_site*1.0d0 / PWDINT%POWER  
        ! it is fine to assume all power generated in fuel assemblies since this is used as a constant scale to quantify power difference
        maxup = 0.0d0
        rmsdif= 0.0d0
        do ii = 1, num_potential_site
            IASS = fuel_pos(ii)
            if(IASS == reload_pos) then
                pwrdif = inverse_avgpwr * (asspwr_in(reload_pos,iload) - refpwr(IASS))
            else
                pwrdif = inverse_avgpwr * (eocpwr(IASS)*scale - refpwr(IASS))
            end if
            if(pwrdif > maxup) then 
                maxup = pwrdif
                maxpos= IASS
            end if
            rmsdif = rmsdif + pwrdif*pwrdif
        end do
        rmsdif = sqrt(rmsdif/num_potential_site)
        
    End Subroutine
end subroutine
