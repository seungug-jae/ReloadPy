Subroutine reload_null(NOUT)
    use m_reload
    use m_utility, only: system_command, insertSort
    implicit none
    integer :: NOUT
    ! local
    integer :: IASS, IX, IY, iload, ipos, num_candidate, num_scratch, ii
    integer, pointer :: candidate_pool(:), candidate_pos(:), scratch_pool(:), scratch_pool2(:,:)
    real(8), pointer :: keffs(:), peaks(:), peakden(:)
    real(8) :: fdum, rmsdif, avgpwr, digit_shift
    character(len=8) :: cdum
    logical :: yes

    if(candidate_range_size > 0) then
        write(NOUT,'(/,"Initital candidate position range:",/,1000I5)') candidate_range_pos
        num_candidate = candidate_range_size
        candidate_pool => candidate_range_pos
    else
        write(NOUT,'(/,"Initital candidate position range: all fuel assemblies")')
        num_candidate = num_potential_site
        candidate_pool => fuel_pos
    end if
    allocate(candidate_pos(num_candidate), keffs(num_candidate), peaks(num_candidate), scratch_pool(num_candidate), &
    &   scratch_pool2(2,num_candidate), peakden(num_candidate))
    dif3d_check_all = .true.  ! force to check all
    
    do iload = 1, num_reload_comp
        ! peak assembly power
        write(NOUT,'(/,20X,"*** Estimating assembly powers (Watt) with reloaded composition No.",I4,1X,A," ***",/)') &
        &   iload, RFUEL(iload)%name
        write(NOUT,'("   Candidate position     Assembly powers from REBUS        Estimated assembly powers after fuel reloading")')
        if(num_hex_sector > 0) then
            write(NOUT,'("    IASS  IRING  IPOS         BOC            EOC               Reloaded         Peak       Peak Assembly",/)')
        else
            write(NOUT,'("    IASS    IX    IY          BOC            EOC               Reloaded         Peak       Peak Assembly",/)')
        end if
        do ipos = 1, num_candidate
            IASS = candidate_pool(ipos)
            if(num_hex_sector > 0) then
                write(NOUT,'(I8,2I6,2X,2(2X,ES13.5),3X,2(4X,ES13.5),3X,I8)') IASS, hexagon_rp(1:2,IASS), &
                    bocpwr(IASS), eocpwr(IASS), asspwr(IASS, iload), pekpwr(IASS,iload), pekass(IASS,iload)
            else
                call reload_ass2xy(IASS, num_xnode, IX, IY)
                write(NOUT,'(I8,2I6,2X,2(2X,ES13.5),3X,2(4X,ES13.5),3X,I8)') IASS, IX, IY, bocpwr(IASS), &
                    eocpwr(IASS), asspwr(IASS, iload), pekpwr(IASS,iload), pekass(IASS,iload)
            end if
        end do

        ! peak node power
        write(NOUT,'(/,20X,"*** Estimating peak power density (W/cc) with reloaded composition No.",I4,1X,A," ***",/)') &
        &   iload, RFUEL(iload)%name
        write(NOUT,'("   Candidate position    Peak power density from DIF3D     Estimated peak power density with reloaded fuel")')
        if(num_hex_sector > 0) then
            write(NOUT,'("    IASS  IRING  IPOS         BOC            EOC             Peak Density    Peak Assembly    Peak Node",/)')
        else
            write(NOUT,'("    IASS    IX    IY          BOC            EOC             Peak Density    Peak Assembly    Peak Node",/)')
        end if
        do ipos = 1, num_candidate
            IASS = candidate_pool(ipos)
            if(num_hex_sector > 0) then
                write(NOUT,'(I8,2I6,2X,2(2X,ES13.5),7X,ES13.5,2I12)') IASS, hexagon_rp(1:2,IASS), &
                    bocden(IASS), eocden(IASS), pekden(IASS, iload), pekpos(:,IASS,iload)
            else
                call reload_ass2xy(IASS, num_xnode, IX, IY)
                write(NOUT,'(I8,2I6,2X,2(2X,ES13.5),7X,ES13.5,2I12)') IASS, IX, IY, &
                    bocden(IASS), eocden(IASS), pekden(IASS, iload), pekpos(:,IASS,iload)
            end if
        end do

        ! calculate assembly worth change
        call reload_reactivity(NOUT, num_candidate, candidate_pool)
        write(NOUT,'(/,20X,"*** Estimating reactivity change due to reloaded composition No.",I4,1X,A," ***")') &
        &   iload, RFUEL(iload)%name 
        write(NOUT,'(" Note: assembly total worth is estimated as k-effective change with composition replacement with ",A8,/)') &
        &   RFUEL(num_reload_comp+1)%name
        if(num_hex_sector > 0) then
            write(NOUT,'("Assembly                  Assembly worth at this position       Added worth",/, &
            &            " index     IRING  IPOS    Burned assembly   Fresh assembly      by fueling",/)')
        else
            write(NOUT,'("Assembly                  Assembly worth at this position       Added worth",/, &
            &            " index       IX    IY     Burned assembly   Fresh assembly      by fueling",/)')
        end if
        do ipos = 1, num_candidate
            IASS = candidate_pool(ipos)
            ! compute worth of burned assembly
            fdum = asswth(IASS,iload) - deltak(IASS,iload)      ! asswth is worth of fresh fuel assembly if reloaded, deltak is added worth
            if(num_hex_sector > 0) then
                write(NOUT,'(I6,3X,2I6,4X,F13.5,5X,F13.5,4X,F13.5)') IASS, hexagon_rp(:,IASS),fdum,asswth(IASS,iload),deltak(IASS,iload)
            else
                call reload_ass2xy(IASS, num_xnode, IX, IY)
                write(NOUT,'(I6,3X,2I6,4X,F13.5,5X,F13.5,4X,F13.5)') IASS, IX, IY, fdum, asswth(IASS,iload), deltak(IASS,iload)
            end if
        end do

        ! perform dif3d check
        if(dif3d_check) then
            ! set up dif3d run environment
            call system_command(NOUT, 'mkdir '//dif3d_dir)  ! make a separate working directory for DIF3D
            inquire(file='NHFLUX', exist=yes)
            if(yes) then
                call reload_dif3d_input(NOUT,'reload.dif3d',1)  ! make dif3d input using interface files
                call system_command(NOUT, 'cp NHFLUX '//dif3d_dir//'/NHFLUX')
            else
                call reload_dif3d_input(NOUT,'reload.dif3d',0)
            end if
            call system_command(NOUT, 'cp reload.dif3d '//dif3d_dir//'/dif3d.inp')
            candidate_pos(1:num_candidate) = candidate_pool(1:num_candidate)
            if(objective_func == 'assembly_power') then
                do ii = 1, num_candidate
                    peaks(ii) = pekpwr(candidate_pos(ii),iload)
                end do
                call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend', eps_power)
            else if(objective_func == 'peak_power') then
                do ii = 1, num_candidate
                    peaks(ii) = pekden(candidate_pos(ii),iload)
                end do
                call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend', eps_power)
            else if(objective_func == 'burnup') then
                do ii = 1, num_candidate
                    peaks(ii) = assbrn(candidate_pos(ii))
                end do
                call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'descend', eps_power)
            else if(objective_func == 'reactivity' .and. use_persent) then
                do ii = 1, num_candidate
                    peaks(ii) = deltak(candidate_pos(ii),iload)
                end do
                call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'descend',-eps_keff)
            else if(objective_func == 'reactivity_slow' .and. use_persent) then
                do ii = 1, num_candidate
                    peaks(ii) = deltak(candidate_pos(ii),iload) - (pekfac(candidate_pos(ii),iload) - eocpekfac)*wgt_pekfac*conv_peak_rho
                end do
                call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'descend',-eps_keff)
            else if(objective_func == 'power_shape') then
                avgpwr = PWDINT%POWER / num_potential_site  ! approximate average fuel assembly power
                digit_shift = 10**(-floor(log10(eps_power_diff)))      ! reasonable assumption: eps_power_diff < 1.0
                do ii = 1, num_candidate
                    IASS = candidate_pos(ii)
                    call reload_power_change(iload, IASS, fdum, ipos, rmsdif)
                    ! since many candidates may have similar rmsdif, sorting first by increasing rmsdif and then by increasing peak assembly power
                    fdum = pekpwr(IASS,iload) / avgpwr      ! should be well below 10.0
                    peaks(ipos) = floor(rmsdif*digit_shift+0.5d0)*10 + fdum
                end do
                call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend', -eps_power)
            else if(objective_func == 'power_deep') then
                do ii = 1, num_candidate
                    IASS = candidate_pos(ii)
                    peaks(ii) = eocpwr(IASS) / refpwr(IASS) - 1.0d0 ! power change relative to reference condition
                end do
                call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend',-eps_power)
            end if
            write(NOUT,'(/," dif3d check will be performed for candidates in the following order:",/,1000I5)') &
            &   candidate_pos(1:num_candidate)
            ! run dif3d with replaced assembly at candidate positions
            call reload_rundif3d(NOUT, num_candidate, candidate_pos, iload, RFUEL(iload),keffs(1:num_candidate),peaks(1:num_candidate),&
            &   scratch_pool(1:num_candidate), peakden(1:num_candidate), scratch_pool2(:,1:num_candidate), num_scratch)    ! now use scratch_pool to record peak assembly position for each reloading
            !call system_command(NOUT, 'rm -r '//dif3d_dir)      ! remove this if like to see DIF3D output

            write(NOUT,'(/,20X,"*** DIF3D calculations with reloaded assembly at candidate positions ***")')
            if(num_hex_sector > 0) then
                write(NOUT,'(/,36X,"Peak assembly integrated power     Peak node-averaged power density    Valid",/,&
            &"Index  IRING  IPOS    Delta_Keff    Value (W)   Assembly Location      Value (W/cc)    Assembly    Node    Candidate?",/)')
            else
                write(NOUT,'(/,36X,"Peak assembly integrated power     Peak node-averaged power density    Valid",/, &
            &"Index    IX    IY     Delta_Keff    Value (W)   Assembly Location      Value (W/cc)    Assembly    Node    Candidate?",/)')
            end if
            do ipos = 1, num_scratch
                IASS = candidate_pos(ipos)
                if(IASS > 0) then
                    cdum = 'Yes'
                else
                    cdum = 'No'
                    IASS = -IASS
                end if
                if(num_hex_sector > 0) then
                    write(NOUT,'(I5,2I6,3X,F11.5,3X,ES13.5,3X,I6,13X,ES13.5,2I10,5X,A)') IASS,hexagon_rp(1:2,IASS),keffs(ipos)-eoc_keff,&
                    &   peaks(ipos),scratch_pool(ipos), peakden(ipos), scratch_pool2(:,ipos), cdum
                else
                    call reload_ass2xy(IASS, num_xnode, IX, IY)
                    write(NOUT,'(I5,2I6,3X,F11.5,3X,ES13.5,3X,I6,13X,ES13.5,2I10,5X,A)') IASS, IX, IY, keffs(ipos)-eoc_keff, &
                    &   peaks(ipos),scratch_pool(ipos), peakden(ipos), scratch_pool2(:,ipos), cdum
                end if
            end do
        end if
    end do
    call system_command(NOUT,'mv  user.NAFLUX  init.NAFLUX')   ! saved as initial guess for next cycle 
    call system_command(NOUT,'rm  user.NHFLUX')   ! will have new NHFLUX in next cycle 

End Subroutine