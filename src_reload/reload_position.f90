Subroutine reload_position(NOUT)
    use m_reload
    use m_utility
    implicit none
    integer :: NOUT
    ! local
    integer :: IASS, IZ, IX, IY, iload, nreg, ireg, MRXY(GEODST%NINTK)
    integer :: reload_pos, reload_pekass, reload_pekpos(2)
    integer :: candidate_pos(num_potential_site), ipos, num_candidate, num_pool, initial_pool_size
    integer :: scratch_pool(num_potential_site), scratch_pool2(2,num_potential_site), num_scratch, ii, kk
    integer,pointer :: initial_pool(:)
    real(8) :: peak_pwrden, peak_asspwr, peak_burnup, peak_keff, peak_dist  ! not necessarily peak value, but quantities corresponding to final reloading position
    character(len=8) :: reload_reg(top_active_node), cdum
    real(8) :: keffs(num_potential_site)   ! keff with reloaded assembly at candidate positions
    real(8) :: peaks(num_potential_site)   ! used as scratch working array at first and then store peak assembly powers obtained with DIF3D calculations
    real(8) :: peakden(num_potential_site) ! peak power density obtained with DIF3D calculation for each candidate reloading position
    real(8) :: maxup(num_potential_site)   ! maximum incurease in assembly power of reloaded core from reference condition
    real(8) :: rmsdif(num_potential_site)  ! root mean square change in assembly power distribution
    logical :: search_failed, update_position, done_with_persent, yes
    real(8) :: fdum, fdum2, tic, toc, avgpwr, digit_shift
    integer :: maxpos
    integer :: iSearch(num_reload_comp)    ! search results for each reload composition, 0 means search failed

    100 format('[reload_position]...',A,1X,I6)
    200 format('[reload_position]...Error: no refueling position satisfies constraint of',1X,A,1X,/,'Program terminated')

    ! power peaking constraint
    if(fixed_peak_power) then  ! using a fixed peak power limit
        if(peak_power_limit < 0.0) peak_power_limit = boc_peak_power_density * (1.0 + peak_tolerance)     ! peak_tolerance initialized in m_reload
        if(assembly_power_limit < 0.0) assembly_power_limit = boc_peak_assembly_power * (1.0 + peak_tolerance)
    else
        peak_power_limit = maxval(eocden)
        assembly_power_limit = maxval(eocpwr)
    end if

    if(candidate_range_size > 0) then
        write(NOUT,'(/,"Initital candidate position range:",/,1000I5)') candidate_range_pos
        initial_pool => candidate_range_pos
        initial_pool_size = candidate_range_size
    else 
        write(NOUT,'(/,"Initital candidate position range: all fuel assemblies")')
        initial_pool => fuel_pos
        initial_pool_size = num_potential_site
    end if

    ! 
    
    ! determine reloading position
    do iload = 1, num_reload_comp
        search_failed = .false.
        write(NOUT,'(/,20X,"*** Estimating assembly powers (Watt) with reloaded composition No.",I4,1X,A," ***",/)') &
        &   iload, RFUEL(iload)%name
        write(NOUT,'("   Candidate position     Assembly powers from REBUS        Estimated assembly powers after fuel reloading")')
        if(num_hex_sector > 0) then
            write(NOUT,'("    IASS  IRING  IPOS         BOC            EOC               Reloaded         Peak       Peak Assembly",/)')
        else
            write(NOUT,'("    IASS    IX    IY          BOC            EOC               Reloaded         Peak       Peak Assembly",/)')
        end if
        do ipos = 1, initial_pool_size
            IASS = initial_pool(ipos)
            if(num_hex_sector > 0) then
                write(NOUT,'(I8,2I6,2X,2(2X,ES13.5),3X,2(4X,ES13.5),3X,I8)') IASS, hexagon_rp(1:2,IASS), &
                    bocpwr(IASS), eocpwr(IASS), asspwr(IASS, iload), pekpwr(IASS,iload), pekass(IASS,iload)
            else
                call reload_ass2xy(IASS, num_xnode, IX, IY)
                write(NOUT,'(I8,2I6,2X,2(2X,ES13.5),3X,2(4X,ES13.5),3X,I8)') IASS, IX, IY, bocpwr(IASS), &
                    eocpwr(IASS), asspwr(IASS, iload), pekpwr(IASS,iload), pekass(IASS,iload)
            end if
        end do

        write(NOUT,'(/,20X,"*** Estimating peak power density (W/cc) with reloaded composition No.",I4,1X,A," ***",/)') &
        &   iload, RFUEL(iload)%name
        write(NOUT,'("   Candidate position    Peak power density from DIF3D     Estimated peak power density with reloaded fuel")')
        if(num_hex_sector > 0) then
            write(NOUT,'("    IASS  IRING  IPOS         BOC            EOC             Peak Density    Peak Assembly    Peak Node",/)')
        else
            write(NOUT,'("    IASS    IX    IY          BOC            EOC             Peak Density    Peak Assembly    Peak Node",/)')
        end if
        do ipos = 1, initial_pool_size
            IASS = initial_pool(ipos)
            if(num_hex_sector > 0) then
                write(NOUT,'(I8,2I6,2X,2(2X,ES13.5),7X,ES13.5,2I12)') IASS, hexagon_rp(1:2,IASS), &
                    bocden(IASS), eocden(IASS), pekden(IASS, iload), pekpos(:,IASS,iload)
            else
                call reload_ass2xy(IASS, num_xnode, IX, IY)
                write(NOUT,'(I8,2I6,2X,2(2X,ES13.5),7X,ES13.5,2I12)') IASS, IX, IY, &
                    bocden(IASS), eocden(IASS), pekden(IASS, iload), pekpos(:,IASS,iload)
            end if
        end do
        ! ================== apply all constraints to narrow down candidate pool
        ! screen with peak power constraint
        write(NOUT,'(/," limit of peak power density (W/cc) =",ES13.5)') peak_power_limit
        write(NOUT,'(" limit of assembly-integrated power (W) =",ES13.5)') assembly_power_limit
        if(manual_reload) then
            num_candidate = 1
            candidate_pos(1) = manual_reload_pos(1) + num_xnode * (manual_reload_pos(2)-1)
            IASS = candidate_pos(1)
            fdum = pekden(IASS,iload) / peak_power_limit - 1.0D0
            fdum2 = pekpwr(IASS,iload) / assembly_power_limit - 1.0D0
            if(fdum > eps_power .or. fdum2 > eps_power) then
                write(NOUT,100) 'Note: peak power limit constraint violated by user selected refueling position',IASS
            else
                write(NOUT,100) 'Note: peak power limit constraint satisfied by user selected refueling position',IASS
            end if
        else
            num_candidate = 0
            do ipos = 1, initial_pool_size
                IASS = initial_pool(ipos)
                fdum = pekden(IASS,iload) / peak_power_limit - 1.0D0
                fdum2 = pekpwr(IASS,iload) / assembly_power_limit - 1.0D0
                if(fdum <= eps_power .and. fdum2 <= eps_power) then
                    num_candidate = num_candidate + 1
                    candidate_pos(num_candidate) = IASS
                end if
            end do
            write(NOUT,'(" reloading positions yielding less or equal peak powers:",/,1000I5)') candidate_pos(1:num_candidate)
            if(num_candidate <= 0) then
                search_failed = .true.
                go to 900
            end if

            ! screen out positions inducing small power increase with refueling (indicating small reactivity addition)
            scratch_pool = candidate_pos
            do ipos = 1, num_candidate
                IASS = candidate_pos(ipos)
                fdum = asspwr(IASS,iload) - eocpwr(IASS)
                ! check on this, 150.0 for 0.5 days, 210.0 for 1.0 day
                if(fdum <= 210.0) candidate_pos(ipos) = 0   ! tentative threshold, should be determined from history later
            end do
            num_pool = num_candidate
            call reduce_candidate_pool(candidate_pos, num_pool, num_candidate)
            write(NOUT,'(" reloading positions with noticeable power increase after refueling:",/,1000I5)') &
            &   candidate_pos(1:num_candidate)
            if(num_candidate <= 0) then     ! roll back to original candidates and do detailed check
                write(NOUT,'(" no candidate passes this check, will skip this screening")')
                candidate_pos = scratch_pool
                num_candidate = num_pool
            end if
        end if

        ! apply discharge burnup constraint
        if(consider_burnup) then
            write(NOUT,'(/," required minimum burnup (MWD/MT) of discharged assembly =",ES13.5)') min_burnup
            if(objective_func == 'burnup' .and. min_burnup > maxval(assbrn)) then
                write(NOUT,100) 'Note: user specified minimum discharge burnup > observed maximum burnup'
                write(NOUT,100) '  since the prior objective is to maximize discharge burnup anyway, this constraint is released'
            else
                do ipos = 1, num_candidate
                    IASS = candidate_pos(ipos)
                    fdum = assbrn(IASS) / min_burnup - 1.0D0
                    if(fdum < -eps_power) candidate_pos(ipos) = 0
                end do
                num_pool = num_candidate
                call reduce_candidate_pool(candidate_pos, num_pool, num_candidate)
                write(NOUT,'(" reloading positions with greater or equal discharge burnups:",/,1000I5)') &
                &   candidate_pos(1:num_candidate)
                if(num_candidate <= 0) then
                    search_failed = .true.
                    go to 900
                end if
            end if
        end if

        ! apply constraints on reloading distance
        if(disperse_reload) then
            write(NOUT,'(" Assembly positions recently refueled:",1000I5)') previous_pos
            if(disperse_mindist < 0.0) then
                write(NOUT,'(" will exclude these from candidate refueling positions")')
                do ipos = 1, num_candidate
                    IASS = candidate_pos(ipos)
                    do kk = 1, num_previous_site
                        if(IASS == previous_pos(kk)) then
                            candidate_pos(ipos) = 0
                            exit
                        end if
                    end do
                end do
            else
                write(NOUT,'(" required minimum distance (* assembly pitch) from the above assemblies =",ES13.5)') disperse_mindist
                allocate(disperse_dist(num_xynod))
                disperse_dist(:) = -1.0d0 
                ! compute minimum distance of each refueling position to previous refueling positions
                call reload_distance(NOUT, num_candidate, candidate_pos)
                ! apply constraint
                scratch_pool = candidate_pos
                do ipos = 1, num_candidate
                    IASS = candidate_pos(ipos)
                    if(disperse_dist(IASS)/disperse_mindist - 1.0d0 < -eps_dist) candidate_pos(ipos) = 0
                end do
            end if
            num_pool = num_candidate
            call reduce_candidate_pool(candidate_pos, num_pool, num_candidate)
            if(num_candidate == 0 .and. disperse_mindist > 0.0) then
                write(NOUT,'(" Since no refueling position has distance from previous refueling position bigger than ",ES13.5, &
                &   " * assembly pitch, this constraint is released to proceed.")') disperse_mindist
                num_candidate = num_pool
                candidate_pos = scratch_pool
            else if(num_candidate > 0) then
                write(NOUT,'(" reloading positions having equal or larger distance:",/,1000I5)') candidate_pos(1:num_candidate)
            else
                write(NOUT,'(" no candidate refueling position left after screening")') 
                search_failed = .true. 
                go to 900
            end if
        end if
        
        ! apply minimum cycle length constraint
        eps_keoc = eps_keff
        
        if(adjust_cycle_length.and..not.fixed_sequence) then    ! note that eps_keoc forced to be consistent with eps_cycle as there is no cycle length iteration here
            min_deltak = (boc_keff - eoc_keff) / cycle_length * min_cycle_length
            !eps_keoc = min_deltak * eps_cycle
            min_deltak = max(min_deltak, 1.d-6)
            write(NOUT,'("Note: required minimum keff increase determined by minimum cycle length =",F13.5)') min_deltak
        end if
        
        ! apply reactivity constraint
        if(consider_keff) then
            if (fixed_sequence) then
                min_deltak = -1._8
                max_deltak = 100000
                mindk = -1.0
                write(NOUT,'("The sequence is fixed by user",F13.5)')
                write(NOUT,'("Note: required minimum keff increase set by =",F13.5)') min_deltak
                write(NOUT,'("Note: required maximum keff increase set by =",F13.5)') max_deltak
            else
                ! if cycle length is not adjusted and min_deltak not specified, set it to 80% of reactivity loss of one depletion cycle
                if(min_deltak < 0.0) min_deltak = 0.8*(boc_keff-eoc_keff)
                if(min_deltak - max_deltak > eps_keff) then 
                    if(adjust_cycle_length) then
                        write(NOUT,'("Error: min_deltak determined by required minimum cycle length =",F13.5)') min_deltak
                    else
                        write(NOUT,'("Error: min_deltak set to 80% x (BOC_keff - EOC_keff) by default =",F13.5)') min_deltak
                    end if
                    write(NOUT,'("  For this problem, it is greater than specified max_deltak =",F13.5, &
                    &   ", a manual setting of min_deltak/max_deltak is required in this case")') max_deltak
                    call abort
                end if
                if(objective_func == 'reactivity' .and. (.not. adjust_cycle_length)) then
                    write(NOUT,100) 'Note: since the prior objective is to maximize reactivity increase anyway, &
                    &the constraint on minimum increase of k-effective by refueling is released'
                    mindk = -1.0
                else
                    mindk = min_deltak
                end if
            end if
            if(use_persent) then
                select case (trim(objective_func))
                case ('assembly_power') ! sort candidate positions in order of increasing peak assembly power after refueling for batch wise run
                    do ipos = 1, num_candidate
                        IASS = candidate_pos(ipos)
                        peaks(ipos) = pekpwr(IASS,iload) 
                    end do
                    call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend', eps_power)
                case ('peak_power') ! sort candidate positions in order of increasing peak power density after refueling for batch wise run
                    do ipos = 1, num_candidate
                        IASS = candidate_pos(ipos)
                        peaks(ipos) = pekden(IASS,iload) 
                    end do
                    call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend', eps_power)
                case ('burnup') ! sort candidate positions in order of decreasing discharge burnup for batch wise run
                    do ipos = 1, num_candidate
                        IASS = candidate_pos(ipos)
                        peaks(ipos) = assbrn(IASS)
                    end do
                    call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'descend', eps_power)
                case ('reactivity') ! sort candidate positions in order of decreasing relative change of assembly powers due to refueling for batch wise run
                    do ipos = 1, num_candidate
                        IASS = candidate_pos(ipos)
                        peaks(ipos) = asspwr(IASS,iload) / eocpwr(IASS) - 1.0d0     ! relative power change due to refueling, indicator of added reactivity
                        !peaks(ipos) = peaks(ipos) * adjoint_flux        ! estimate reactivity change 
                    end do
                    call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'descend', -eps_power)
                case ('power_shape') ! sort candidate positions in order of increasing deviation of power distribution
                    write(NOUT,'(/,20X,"*** Estimating change in assembly power of reloaded core compared to reference condition ***")')
                    write(NOUT,'(20X,"    (power change relative to the core-averaged fuel assembly power is given below)",/)')
                    write(NOUT,'("candidate reload position   RMS change       Max. increase    Assembly seeing    Corresponding peak   Sorting")')
                    write(NOUT,'("   IASS    IRING   IPOS     /Avg power (%)   /Avg power (%)   maximum increase   assembly power (W)   metric ")')
                    avgpwr = PWDINT%POWER / num_potential_site  ! approximate average fuel assembly power
                    digit_shift = 10**(-floor(log10(eps_power_diff)))      ! reasonable assumption: eps_power_diff < 1.0
                    do ipos = 1, num_candidate
                        IASS = candidate_pos(ipos)
                        call reload_power_change(iload, IASS, maxup(ipos), maxpos, rmsdif(ipos))
                        ! since many candidates may have similar rmsdif, sorting first by increasing rmsdif and then by increasing peak assembly power
                        fdum = pekpwr(IASS,iload) / avgpwr      ! should be well below 10.0
                        peaks(ipos) = floor(rmsdif(ipos)*digit_shift+0.5d0)*10 + fdum
                        ! edit
                        write(NOUT,'(I7,2X,2I7,5X,2(F10.2,7X),I10,9X,ES13.5,6X,F10.3)') IASS, hexagon_rp(:,IASS), rmsdif(ipos)*100, &
                            maxup(ipos)*100, maxpos, pekpwr(IASS,iload), peaks(ipos)
                    end do
                    call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend', -eps_power)
                case ('power_deep')  ! sort candidate positions in order of increasing relative assembly power (Power_eoc/Power_ref) for batch wise run
                    write(NOUT,'(/,20X,"*** Relative change in assembly power at EOC compared to reference condition ***",/)')
                    write(NOUT,'("candidate     IRING   IPOS     ref. power    EOC power     Relative change (%)")')
                    do ipos = 1, num_candidate
                        IASS = candidate_pos(ipos)
                        peaks(ipos) = eocpwr(IASS) / refpwr(IASS) - 1.0d0 ! power change relative to reference power
                        write(NOUT,'(I7,7X,I4,4X,I4,2ES14.5,5X,F10.2)') &
                        &   IASS, hexagon_rp(:,IASS), refpwr(IASS), eocpwr(IASS), peaks(ipos)*100
                    end do
                    call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend', -eps_power_diff)
                case default
                    write(NOUT,100) 'Note: batch wise PERSENT run not implemented for this objective_func,&
                    & will estimate reactivity changes in ONE PERSENT run'
                    num_persent_case = num_candidate
                end select

                write(NOUT,'(/," perturbation calculations will be performed for candidates in the following order:",/,1000I5)') &
                &   candidate_pos(1:num_candidate)
                ! print table header
                write(NOUT,'(/,20X,"*** Estimating reactivity change due to reloaded composition No.",I4,1X,A," ***")') &
                &   iload, RFUEL(iload)%name
                if(consider_worth) then
                    write(NOUT,'(" Note: assembly worth is estimated as k-effective change between two composition replacement with ",&
                    &   A8," and ",A8,/)') RFUEL(iload)%name, RFUEL(num_reload_comp+1)%name
                else
                    write(NOUT,'(" Note: assembly worth is not calculated (only net change in keff is determined)",/)')
                end if
                if(num_hex_sector > 0) then
                    write(NOUT,'("Assembly                  Assembly worth at this position       Added worth",/, &
                    &            " index     IRING  IPOS    Burned assembly   Fresh assembly      by fueling",/)')
                else
                    write(NOUT,'("Assembly                  Assembly worth at this position       Added worth",/, &
                    &            " index       IX    IY     Burned assembly   Fresh assembly      by fueling",/)')
                end if
                ! work with PERSENT
                num_scratch = num_candidate
                scratch_pool(1:num_scratch)= candidate_pos(1:num_candidate)
                done_with_persent = .false.
                
                persent_time(1:3) = 0.0d0
                do while (.not. done_with_persent)
                    if(num_persent_case > 0 .and. num_persent_case < num_scratch) then
                        ! shrink candidate pool and save untested cases
                        num_candidate = num_scratch - num_persent_case
                        candidate_pos(1:num_persent_case) = scratch_pool(1:num_persent_case)
                        scratch_pool(1:num_candidate) = scratch_pool(num_persent_case+1:num_scratch)
                        scratch_pool(num_candidate+1:num_scratch) = 0
                        num_scratch = num_candidate
                        num_candidate = num_persent_case
                    else
                        num_candidate = num_scratch
                        candidate_pos(1:num_candidate) = scratch_pool(1:num_candidate)
                        num_scratch = 0
                        done_with_persent = .true.
                    end if

                    ! check reactivity worth of selected candidates
                    call CPU_TIME(tic)
                    call reload_reactivity(NOUT, num_candidate, candidate_pos) 
                    call CPU_TIME(toc)
                    persent_time(1) = persent_time(1) + toc - tic                   
                    do ipos = 1, num_candidate
                        IASS = candidate_pos(ipos)
                        if(consider_worth) then
                            fdum = asswth(IASS,iload) - deltak(IASS,iload)      ! asswth is worth of fresh fuel assembly if reloaded, deltak is added worth
                            fdum2= asswth(IASS,iload)
                        else
                            fdum = 0.0d0; fdum2 = 0.0d0
                        end if
                        if(num_hex_sector > 0) then
                            write(NOUT,'(I6,3X,2I6,4X,F13.5,5X,F13.5,4X,F13.5)') IASS, hexagon_rp(:,IASS),fdum,fdum2,deltak(IASS,iload)
                        else
                            call reload_ass2xy(IASS, num_xnode, IX, IY)
                            write(NOUT,'(I6,3X,2I6,4X,F13.5,5X,F13.5,4X,F13.5)') IASS, IX, IY, fdum,fdum2,deltak(IASS,iload)
                        end if
                        ! apply reactivity constraints
                        if(mindk > 0.0 .and. deltak(IASS,iload) - mindk < -eps_keoc) candidate_pos(ipos) = 0  ! eps_keff set in m_reload
                        if(deltak(IASS,iload) - max_deltak > eps_keff) candidate_pos(ipos) = 0
                        if(consider_worth .and. asswth(IASS,iload) - max_worth > eps_keff) candidate_pos(ipos) = 0
                    end do
                    num_pool = num_candidate
                    call reduce_candidate_pool(candidate_pos, num_pool, num_candidate)
                    if(num_candidate > 0) then
                        done_with_persent = .true.
                    else if(num_candidate == 0 .and. num_scratch == 0) then
                        search_failed = .true.
                        done_with_persent = .true.
                    end if
                end do
                call system_command(NOUT,'mv  user.NAFLUX  init.NAFLUX')   ! saved as initial guess for next cycle 
                call system_command(NOUT,'rm  user.NHFLUX')   ! will have new NHFLUX in next cycle 
                ! print PERSENT timing (system command line execution time cannot be included)
                !write(NOUT,'(/,"Time cost (sec) of PERSENT run:",   &
                !&            /,"    total cost          =",ES13.5,  &
                !&            /,"    calculation         =",ES13.5,  &
                !&            /,"    reading output      =",ES13.5)') persent_time(1:3)
            end if
            
            write(NOUT,'(/," K-effective at BOC =",F13.5, &
            &    /," K-effective at EOC =",F13.5, &
            &    /," cycle length (days) =",ES13.5, &
            &    /," allowed maximum increase of k-effective after refueling =",F13.5, &
            &    /," required minimum increase of k-effective by refueling   =",F13.5)') &
            &    boc_keff, eoc_keff, cycle_length, max_deltak, min_deltak
            if(consider_worth) write(NOUT,'(" allowed maximum reactivity worth of reloaded assembly   =",F13.5)') max_worth
            if(use_persent) then
                write(NOUT,'(" reloading positions meeting these reactivity constraints:",/,1000I5)') candidate_pos(1:num_candidate)
                if(search_failed) go to 900
            else
                write(NOUT,'(" Since PERSENT is not enabled, the constraints of reactivity will be applied in DIF3D check")')
            end if
        end if

        ! ====================== determine final position of refueling based on objective_func
        if(.not. dif3d_check) then
            reload_pos = candidate_pos(1)
            select case (trim(objective_func))
            case ('assembly_power')
                peak_asspwr= pekpwr(reload_pos,iload)
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
                    fdum = pekpwr(IASS,iload) / peak_asspwr - 1.0D0
                    if(fdum < -eps_power) then
                        reload_pos = IASS
                        peak_asspwr = pekpwr(IASS,iload)
                    end if
                end do
            case ('peak_power')
                peak_pwrden = pekden(reload_pos,iload)
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
                    fdum = pekden(IASS,iload) / peak_pwrden - 1.0D0
                    if(fdum < -eps_power) then
                        reload_pos = IASS
                        peak_pwrden = pekden(IASS,iload)
                    end if
                end do
            case ('burnup')
                peak_burnup = assbrn(reload_pos)
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
                    fdum = assbrn(IASS) / peak_burnup - 1.0D0
                    if(fdum > eps_power) then
                        reload_pos = IASS
                        peak_burnup = assbrn(IASS)
                    end if
                end do
            case ('reactivity')
                peak_keff = deltak(reload_pos,iload)
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
                    if(deltak(IASS,iload) - peak_keff > eps_keff) then
                        reload_pos = IASS
                        peak_keff = deltak(IASS,iload)
                    else if(abs(deltak(IASS,iload) - peak_keff) <= eps_keff) then
                        if(assbrn(IASS)/assbrn(reload_pos) - 1.0d0 > eps_power) then
                            reload_pos = IASS
                            peak_keff = deltak(IASS,iload)
                        end if
                    end if
                end do
            case ('power_shape')
                fdum = rmsdif(1)
                fdum2 = pekpwr(candidate_pos(1),iload)
                update_position = .false.
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
                    if((fdum - rmsdif(ipos)) > eps_power_diff) then 
                        update_position = .true.
                    else if(abs(fdum - rmsdif(ipos)) < eps_power_diff) then
                        if(fdum2 / pekpwr(IASS,iload) - 1.0d0 > eps_power) update_position = .true.
                    end if
                    if(update_position) then
                        reload_pos = IASS
                        fdum = rmsdif(ipos)
                        fdum2 = pekpwr(IASS,iload)
                    end if
                end do
            case ('power_deep')     ! Moltex suggested option
                fdum2 = eocpwr(reload_pos) / refpwr(reload_pos)
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
                    fdum = eocpwr(IASS) / refpwr(IASS)
                    if(fdum - fdum2 < -eps_power_diff) then
                        reload_pos = IASS
                        fdum2 = fdum
                    end if
                end do
            case default
                write(NOUT,100) 'Error: unrecognized objective_func '//objective_func
                call abort('')
            end select
            peak_pwrden = pekden(reload_pos,iload)
            peak_asspwr = pekpwr(reload_pos,iload)
            reload_pekass = pekass(reload_pos,iload)
            reload_pekpos = pekpos(1:2,reload_pos,iload)
            if(consider_keff .and. use_persent) peak_keff = eoc_keff + deltak(reload_pos,iload)
        end if
        
        ! =================== final check with dif3d calculation
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
            else if(objective_func == 'power_shape') then
                write(NOUT,100) 'Warning: DIF3D check for objective function <power_shape> not completed, &
                &relative power change from reference based on DIF3D results not added'
                call abort 

                do ii = 1, num_candidate
                    peaks(ii) = abs(maxup(ii))
                end do
                call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend', -eps_power)
            else if(objective_func == 'power_deep') then
                do ii = 1, num_candidate
                    IASS = candidate_pos(ii)
                    peaks(ii) = eocpwr(IASS) / refpwr(IASS) - 1.0d0 ! power change relative to reference condition
                end do
                call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend',-eps_power)
            end if

            write(NOUT,'(" dif3d check will be performed for candidates in the following order:",/,1000I5)') &
            &   candidate_pos(1:num_candidate)
            ! run dif3d with replaced assembly at candidate positions (from most to least probable), stop when all constraints are ensured
            call reload_rundif3d(NOUT,num_candidate,candidate_pos,iload,RFUEL(iload),keffs(1:num_candidate),peaks(1:num_candidate),&
            &   scratch_pool(1:num_candidate), peakden(1:num_candidate), scratch_pool2(:,1:num_candidate), num_scratch)    ! now use scratch_pool to record peak assembly position for each reloading
            !call system_command(NOUT, 'rm -r '//dif3d_dir)      ! remove this if like to see DIF3D output

            write(NOUT,'(/,20X,"*** DIF3D calculations with reloaded assembly at candidate positions ***")')
            if(num_hex_sector > 0) then
                write(NOUT,'(/,36X,"Peak assembly integrated power     Peak node-averaged power density    Valid",/,&
            &"Index  IRING  IPOS    Delta_Keff    Value (W)   Assembly Location      Value (W/cc)    Assembly    Node    Candidate?",/)')
            else
                write(NOUT,'(/,36X,"Peak assembly integrated power     Peak node-averaged power density    Valid",/,&
            &"Index    IX    IY     Delta_Keff    Value (W)   Assembly Location      Value (W/cc)    Assembly    Node    Candidate?",/)')
            end if
            num_candidate = num_scratch
            peak_asspwr = peaks(1) * 1.5
            peak_pwrden = peakden(1) * 1.5    ! make sure not miss the first one
            peak_burnup = -1.0
            peak_keff   = keffs(1) * 0.5      ! make sure not miss the first one
            do ipos = 1, num_candidate
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
                if(candidate_pos(ipos) < 0) then
                    candidate_pos(ipos) = 0
                    num_scratch = num_scratch - 1
                    cycle
                end if
                
                if(manual_reload) then
                    update_position = .true.
                else
                    update_position = .false.
                    ! to minimize peak assembly power first (equivalent to minimize radial peaking)
                    if(objective_func == 'assembly_power') then
                        fdum = peaks(ipos) / peak_asspwr - 1.0D0 
                        if(fdum < -eps_power) then
                            update_position = .true.
                        else if(abs(fdum) <= eps_power) then  ! if multiple candidates give same peak
                            fdum = peakden(ipos) / peak_pwrden - 1.0D0
                            if(fdum < -eps_power) then 
                                update_position = .true.
                            else if(abs(fdum) <= eps_power) then 
                                fdum = keffs(ipos) - peak_keff
                                if(fdum > eps_keff) then
                                    update_position = .true.
                                else if (abs(fdum) <= eps_keff .and. consider_burnup) then
                                    fdum = assbrn(IASS) / assbrn(reload_pos) - 1.0D0
                                    if(fdum > eps_power) update_position = .true.
                                end if
                            end if
                        end if
                    end if

                    ! to minimize peak power density first
                    if(objective_func == 'peak_power') then
                        fdum = peakden(ipos) / peak_pwrden - 1.0D0 
                        if(fdum < -eps_power) then
                            update_position = .true.
                        else if(abs(fdum) <= eps_power) then  ! if multiple candidates give same peak
                            fdum = peaks(ipos) / peak_asspwr - 1.0D0
                            if(fdum < -eps_power) then
                                update_position = .true.
                            else if(abs(fdum) <= eps_power) then
                                fdum = keffs(ipos) - peak_keff
                                if(fdum > eps_keff) then
                                    update_position = .true.
                                else if (abs(fdum) <= eps_keff .and. consider_burnup) then
                                    fdum = assbrn(IASS) / assbrn(reload_pos) - 1.0D0
                                    if(fdum > eps_power) update_position = .true.
                                end if
                            end if
                        end if
                    end if

                    ! to maximize burnup
                    if(objective_func == 'burnup') then 
                        if(peak_burnup < 0.0) then
                            peak_burnup = assbrn(IASS)
                            reload_pos = IASS
                            peak_pwrden = peakden(ipos)
                            reload_pekass = scratch_pool(ipos)
                            reload_pekpos = scratch_pool2(:,ipos)
                            peak_keff = keffs(ipos)
                            if(disperse_reload .and. disperse_mindist > 0.0) peak_dist = disperse_dist(IASS)
                        else
                            fdum = assbrn(IASS) / peak_burnup - 1.0d0
                            if(abs(fdum) < eps_power) then  ! multiple candidates have same burnup
                                fdum = keffs(ipos) - peak_keff
                                if(fdum > eps_keff) then
                                    update_position = .true.
                                else if(abs(fdum) <= eps_keff .and. disperse_reload .and. disperse_mindist > 0.0) then
                                    fdum = disperse_dist(IASS) / peak_dist - 1.0
                                    if(fdum > eps_dist) update_position = .true.
                                end if
                                if(.not. update_position) then
                                    fdum = peakden(ipos) / peak_pwrden - 1.0d0
                                    if(fdum < -eps_power) update_position = .true.
                                end if                            
                            end if
                        end if
                    end if

                    ! to maximize reactivity
                    if(objective_func == 'reactivity') then
                        fdum = keffs(ipos) - peak_keff
                        if(fdum > eps_keff) then
                            update_position = .true.
                        else if(abs(fdum) <= eps_keff) then  ! peak candidate with smaller peak
                            fdum = peakden(ipos) / peak_pwrden - 1.0d0
                            if(fdum < -eps_power) then
                                update_position = .true.
                            else if(abs(fdum) <= eps_power) then  ! peak candidate with higher burnup
                                if(consider_burnup) then
                                    fdum = assbrn(IASS) / assbrn(reload_pos) - 1.0d0
                                    if(fdum > eps_power) update_position = .true.
                                end if
                            end if
                        end if
                        if(.not. update_position .and. disperse_reload .and. disperse_mindist > 0.0) then
                            fdum = disperse_dist(IASS) / peak_dist - 1.0d0
                            if(fdum > eps_dist) update_position = .true.
                        end if
                    end if
                end if

                if(update_position) then
                    reload_pos = candidate_pos(ipos)
                    peak_asspwr = peaks(ipos)
                    peak_pwrden = peakden(ipos)
                    reload_pekass = scratch_pool(ipos)
                    reload_pekpos = scratch_pool2(:,ipos)
                    peak_keff = keffs(ipos)
                    if(disperse_reload .and. disperse_mindist > 0.0) peak_dist = disperse_dist(reload_pos)
                end if
            end do
            if(num_scratch == 0) search_failed = .true.
        end if
        
900     continue
        ! summary edits
        write(NOUT,'(/,20X,"*** Summary of reloading search for composition No.",I4,1X,A," ***",/)') iload, RFUEL(iload)%name
        write(NOUT,'(" Constraints considered for reloading: peak power limit, assembly power limit")',advance='NO')
        if(consider_burnup) write(NOUT,'(", discharge burnup")',advance='NO')
        if(consider_keff .or. consider_worth) write(NOUT,'(", reactivity")',advance='NO')
        if(adjust_cycle_length) write(NOUT,'(", cycle length")',advance='NO')
        if(disperse_reload) write(NOUT,'(", assembly distance")',advance='NO')
        if(manual_reload) write(NOUT,'(", user selection")', advance='NO')
        
        if(search_failed) then
            write(NOUT,'(/," search for reloading position: FAILURE, no position meets all constraints",/)')
            write(0,'("reload.x error: search for reloading position failed, no position meets all constraints")')
            iSearch(iload) = 0
        else
            write(NOUT,'(/," search for reloading position: SUCCESS",/)')
            if(objective_func == 'reactivity') then
                write(NOUT,'(" reloading assembly position is chosen to have maximum reactivity increase")')
            else if(objective_func == 'burnup') then
                write(NOUT,'(" reloading assembly position is chosen to have maximum discharge burnup")')
            else if(objective_func == 'power_shape') then
                write(NOUT,'(" reloading assembly position is chosen to have close assembly power distribution as reference")')
            else if(objective_func == 'power_deep') then
                write(NOUT,'(" reloading assembly that has minimum assembly power relative to user specified reference value")')
            else
                write(NOUT,'(" reloading assembly position is chosen to have minimum peak power")')
            end if
            write(NOUT,'(" reloading position (assembly index):",I6)') reload_pos
            write(NOUT,'("   corresponding location (IASS) giving peak assembly-integrated power:",I6)') reload_pekass
            write(NOUT,'("   corresponding peak assembly-integrated power (W):",ES13.5)') peak_asspwr
            write(NOUT,'("   corresponding location (IASS, IZ) giving peak power density:",2I6)') reload_pekpos(1:2)
            write(NOUT,'("   corresponding peak power density (W/cc):",ES13.5)') peak_pwrden
            write(NOUT,'("   discharge burnup (MWD/MT) of reloaded assembly:",ES13.5)') assbrn(reload_pos)
            if(consider_keff) write(NOUT,'("   k-effective with this reloaded assembly:",f10.5)') peak_keff
            if(adjust_cycle_length) then
                ! target_eoc_keff for next cycle could be adjusted if deprate_ratio /= 1.0
                target_eoc_keff = peak_keff - deprate_ratio * (peak_keff - target_eoc_keff)
                fdum = cycle_length * (peak_keff - target_eoc_keff) / (boc_keff - eoc_keff)
                write(NOUT,'("   targeted keff at EOC of next cycle:",F13.5)') target_eoc_keff
                write(NOUT,'("   estimated cycle length (days) for next cycle:",ES13.5)') fdum
            end if
            ! mark regions to be modified in next REBUS calculation
            call reload_ass2xy(reload_pos, num_xnode, IX, IY)
            MRXY(:) = GEODST%MR(IX,IY,:)
            nreg = 0
            do IZ = bot_active_node, top_active_node    ! assumed fuel nodes located at the same axial section in all fuel assemblies
                ireg = MRXY(IZ)
                if(IZ > bot_active_node .and. ireg == MRXY(IZ-1)) cycle
                nreg = nreg + 1
                reload_reg(nreg) = LABELS%REGNAM(ireg)
            end do
            iSearch(iload) = reload_pos
            if (measure_uncertainty) then
                call reload_uncertainty(NOUT, reload_pos)
                call system_command(NOUT,'mv  user.NAFLUX  init.NAFLUX')   ! saved as initial guess for next cycle 
                call system_command(NOUT,'rm  user.NHFLUX')   ! will have new NHFLUX in next cycle 
            end if
        end if
    end do

    ! place holder for now (pick one among multiple reload composition and position)
    write(NOUT,'(/,20X,"*** Summary of reloaded assembly ***",/)')
    write(NOUT,'("reloaded         reload_position     reload_regions",/,"composition      IASS   IX   IY      names")')
    do iload = 1, num_reload_comp
        if(iSearch(iload) > 0) &
        &   write(NOUT,'(I2,3X,A8,3X,3I5,6X,100(A8,1X))') iload, RFUEL(iload)%name, reload_pos, IX, IY, reload_reg(1:nreg)
    end do

End Subroutine

! ----------------------------------------------------------------------
Subroutine reload_ass2xy(IASS, NX, IX, IY)
    integer :: IASS, NX, IX, IY
    IY = (IASS+NX-1) / NX
    IX = IASS - NX * (IY-1)
End Subroutine

! ----------------------------------------------------------------------
Subroutine reduce_candidate_pool(pool, num_pool, num_candidate)
    implicit none 
    integer, intent(in)    :: num_pool
    integer, intent(out)   :: num_candidate
    integer, intent(inout) :: pool(num_pool)
    integer :: ii
    num_candidate = 0
    do ii = 1, num_pool 
        if(pool(ii) > 0) then
            num_candidate = num_candidate + 1
            pool(num_candidate) = pool(ii)
        end if
    end do
    if(num_candidate < num_pool) pool(num_candidate+1:num_pool) = 0
End Subroutine

