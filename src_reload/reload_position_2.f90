module find_position
  implicit none
  
    type candidate_info_type
        logical :: lpassed = .FALSE.
        real(8) :: keffs, peaks, peakden
        real(8) :: refpwr ! reference power used when obj function is 'power_deep'
        integer :: iass ! assembly index
        real(8) :: assbrn ! assembly burnup
        real(8) :: disperse_dist ! disperse distance
        real(8) :: asswth ! assembly reload worth (only after persent calc. / i.e. only for candidate passed for first exam)
        ! data about power shape after reload this position
        real(8) :: maxup, rmsdif
        integer :: maxpos 
        ! assembly powers from REBUS (previous cycle)
        real(8) :: bocpwr, bocden
        real(8) :: eocpwr, eocden
        ! Estimated assembly powers after fuel reloading
        real(8) :: asspwr ! reloaded assembly
        real(8) :: pekpwr ! peak assembly power 
        real(8) :: pekfac ! radial peaking factor (assembly)
        integer :: pekass ! index of peak assembly
        !    peak node
        real(8) :: pekden ! peak node power density
        integer :: pekpos(2)
    end type
    type(candidate_info_type), pointer :: cand_info(:,:)
    integer, pointer :: iass2icand(:)
    integer :: max_iass 
  contains
  
Subroutine reload_position_2(NOUT)
    use m_reload
    use m_utility
    implicit none
    integer :: NOUT
    ! local
    integer :: IASS, IZ, IX, IY, iload, nreg, ireg, MRXY(GEODST%NINTK)
    integer :: reload_pos, reload_pekass, reload_pekpos(2)
    integer :: candidate_pos(num_potential_site), ipos, num_candidate, num_pool, initial_pool_size
    integer :: scratch_pool(num_potential_site), scratch_pool2(2,num_potential_site), num_scratch, ii, kk
  integer :: num_buf_persent, buf_persent(num_potential_site)
    !
    integer,pointer :: initial_pool(:)
    real(8) :: peak_pwrden, peak_asspwr, peak_burnup, peak_keff, peak_dist  ! not necessarily peak value, but quantities corresponding to final reloading position
    real(8) :: peak_pekfac
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
    !
    logical :: lselected ! select the candidate
  logical :: linit_persent
    integer :: nsearch_iter
    INTEGER :: npersent_out_stack, ndif3d_out_stack
    Character*256 :: persent_out_stack(num_potential_site)
    Character*256 :: dif3d_out_stack(num_potential_site)
    ! 
    integer :: passlist_power(num_potential_site), npassed_power
    integer :: passlist_burnup(num_potential_site), npassed_burnup
    integer :: passlist_distance(num_potential_site), npassed_distance
    real(8) :: margin_weight
    logical :: lpass_all
    logical :: lpass_power, lpass_minpower, lpass_burnup, lpass_distance, lpass_reactivity
    logical :: lrelease_burnup
    ! buffer for candidate --------------------------------------------
    integer :: candidate_buf(num_potential_site)
    type(candidate_info_type), pointer :: tempcand
    ! ----------------------------------------------------------------
    100 format('[reload_position]...',A,1X,I6)
    200 format('[reload_position]...Error: no refueling position satisfies constraint of',1X,A,1X,/,'Program terminated')

    if(candidate_range_size > 0) then
        write(NOUT,'(/,"Initital candidate position range:",/,1000I5)') candidate_range_pos
        initial_pool => candidate_range_pos
        initial_pool_size = candidate_range_size
    else 
        write(NOUT,'(/,"Initital candidate position range: all fuel assemblies")')
        initial_pool => fuel_pos
        initial_pool_size = num_potential_site
    end if
    
    max_iass= initial_pool(1)
    do ipos = 2, initial_pool_size
        max_iass = MAX(max_iass, initial_pool(ipos))
    END DO
    if (associated(iass2icand)) deallocate(iass2icand)
    ALLOCATE(iass2icand(max_iass))
    do ipos = 1, initial_pool_size
        iass2icand(initial_pool(ipos)) = ipos
    END DO
    ! power peaking constraint
    if(fixed_peak_power) then  ! using a fixed peak power limit
        if(peak_power_limit < 0.0) peak_power_limit = boc_peak_power_density * (1.0 + peak_tolerance)     ! peak_tolerance initialized in m_reload
        if(assembly_power_limit < 0.0) assembly_power_limit = boc_peak_assembly_power * (1.0 + peak_tolerance)
    else
        peak_power_limit = maxval(eocden)
        assembly_power_limit = maxval(eocpwr)
    end if
    if (associated(cand_info)) deallocate(cand_info)
    allocate(cand_info(num_potential_site, num_reload_comp))
    ! determine reloading position
    do iload = 1, num_reload_comp
        lselected = .FALSE.
        !------------------------------------------------------------------------------------------------------------------------------
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
        end do ! ipos
            
        
        write(NOUT, *) 
        write(NOUT,'(" initial margin for power limit =", ES13.5)') margin_power_limit
        write(NOUT,'(" limit of peak power density (W/cc) =",ES13.5)') peak_power_limit
        write(NOUT,'(" limit of assembly-integrated power (W) =",ES13.5)') assembly_power_limit
        if(consider_burnup) then
            lrelease_burnup = .FALSE.
            write(NOUT,'(" required minimum burnup (MWD/MT) of discharged assembly =",ES13.5)') min_burnup
            if(objective_func == 'burnup' .and. min_burnup > maxval(assbrn)) then
                write(NOUT,100) 'Note: user specified minimum discharge burnup > observed maximum burnup'
                write(NOUT,100) '  since the prior objective is to maximize discharge burnup anyway, this constraint is released'
                lrelease_burnup = .TRUE.
            end if
        end if
        if(disperse_reload) then
            write(NOUT,'(" Assembly positions recently refueled:",1000I5)') previous_pos
            if(disperse_mindist >= 0.0) then
                write(NOUT,'(" required minimum distance (* assembly pitch) from the above assemblies =",ES13.5)') disperse_mindist
                allocate(disperse_dist(num_xynod))
                disperse_dist(:) = -1.0d0 
                ! compute minimum distance of each refueling position to previous refueling positions
                call reload_distance(NOUT, num_candidate, candidate_pos)
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
            
            write(NOUT,'(/," K-effective at BOC =",F13.5, &
            &    /," K-effective at EOC =",F13.5, &
            &    /," cycle length (days) =",ES13.5, &
            &    /," allowed maximum increase of k-effective after refueling =",F13.5, &
            &    /," required minimum increase of k-effective by refueling   =",F13.5)') &
            &    boc_keff, eoc_keff, cycle_length, max_deltak, min_deltak
            if(consider_worth) write(NOUT,'(" allowed maximum reactivity worth of reloaded assembly   =",F13.5)') max_worth
            if(consider_worth) then
                write(NOUT,'(" Note: assembly worth is estimated as k-effective change between two composition replacement with ",&
                &   A8," and ",A8,/)') RFUEL(iload)%name, RFUEL(num_reload_comp+1)%name
            else
                write(NOUT,'(" Note: assembly worth is not calculated (only net change in keff is determined)",/)')
            end if
        end if
        
        ! ================== in case manual reload
        if (manual_reload) then
            ! screen with peak power constraint
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
            ! screen with peak power constraint
            lpass_power = pass_fail(nout, tempcand%pekden, peak_power_limit, eps_power, 1) .AND. pass_fail(nout, tempcand%pekpwr, assembly_power_limit, eps_power, 1)
            lpass_minpower = (tempcand%asspwr - tempcand%eocpwr) >= 210.0
            IF (lpass_power .AND. lpass_minpower) THEN
                write(NOUT,100) 'Note: peak power limit constraint satisfied by user selected refueling position',IASS
            ELSE
                write(NOUT,100) 'Note: peak power limit constraint violated by user selected refueling position',IASS
            END IF
            ! apply discharge burnup constraint
            lpass_burnup = .TRUE.
            if (consider_burnup .AND. .NOT. lrelease_burnup) then
                lpass_burnup = pass_fail(nout, tempcand%assbrn, min_burnup, eps_power, 2)
            end if
            IF (consider_burnup) THEN
                IF (lpass_burnup) THEN
                    write(NOUT,100) 'Note: burnup limit constraint satisfied by user selected refueling position',IASS
                ELSE
                    write(NOUT,100) 'Note: burnup limit constraint violated by user selected refueling position',IASS
                END IF
            END IF
            ! apply disperse distance
            IF (disperse_reload) then
                if(disperse_mindist < 0.0) then
                    lpass_distance = .NOT. lReaptedASM(iass, previous_pos, num_previous_site)
                else
                    lpass_distance = pass_fail(nout, tempcand%disperse_dist, disperse_mindist, eps_dist, 2)
                end if
                IF (lpass_distance) THEN
                    write(NOUT,100) 'Note: disperse distance limit constraint satisfied by user selected refueling position',IASS
                ELSE
                    write(NOUT,100) 'Note: disperse distance limit constraint violated by user selected refueling position',IASS
                END IF
            end if
            
            cand_info(1, iload)%iass   = IASS
            cand_info(1, iload)%bocpwr = bocpwr(IASS)
            cand_info(1, iload)%eocpwr = eocpwr(IASS)
            cand_info(1, iload)%asspwr = asspwr(IASS, iload)
            cand_info(1, iload)%pekpwr = pekpwr(IASS, iload)
            cand_info(1, iload)%pekfac = pekfac(IASS, iload)
            cand_info(1, iload)%pekass = pekass(IASS, iload)
            cand_info(1, iload)%bocden = bocden(IASS)
            cand_info(1, iload)%eocden = eocden(IASS)
            cand_info(1, iload)%pekden = pekden(IASS, iload)
            cand_info(1, iload)%pekpos = pekpos(:,IASS,iload)
            cand_info(1, iload)%assbrn = assbrn(IASS)
            if(objective_func == 'power_deep') cand_info(1, iload)%refpwr = refpwr(IASS)
            if (disperse_reload .AND. disperse_mindist >= 0.0) then
                cand_info(1, iload)%disperse_dist = disperse_dist(iass)
            end if
            iass2icand(iass) = 1
        else        ! ================== apply all constraints to narrow down candidate pool
            !=============================================== first exam (power, burnup, disperse distance)
            ! store the information
            WRITE(NOUT, '(A)') '***** USER EDIT FOR DEVELOPEMENT 1'
            WRITE(NOUT, '(A5, 8A13, A10, 2A13)') 'iass', 'bocpwr', 'eocpwr', 'asspwr', 'pekpwr', 'pekass', &
                                                         'bocden', 'eocden', 'pekden', 'pekpos', 'assbrn', 'rmsdif'
            do ipos = 1, initial_pool_size
                IASS = initial_pool(ipos)
                cand_info(ipos, iload)%iass   = IASS
                cand_info(ipos, iload)%bocpwr = bocpwr(IASS)
                cand_info(ipos, iload)%eocpwr = eocpwr(IASS)
                cand_info(ipos, iload)%asspwr = asspwr(IASS, iload)
                cand_info(ipos, iload)%pekpwr = pekpwr(IASS, iload)
                cand_info(ipos, iload)%pekfac = pekfac(IASS, iload)
                cand_info(ipos, iload)%pekass = pekass(IASS, iload)
                cand_info(ipos, iload)%bocden = bocden(IASS)
                cand_info(ipos, iload)%eocden = eocden(IASS)
                cand_info(ipos, iload)%pekden = pekden(IASS, iload)
                cand_info(ipos, iload)%pekpos = pekpos(:,IASS,iload)
                cand_info(ipos, iload)%assbrn = assbrn(IASS)
                call reload_power_change(iload, IASS, cand_info(ipos, iload)%maxup, cand_info(ipos, iload)%maxpos, cand_info(ipos, iload)%rmsdif)
                if(objective_func == 'power_deep') cand_info(ipos, iload)%refpwr = refpwr(IASS)
                if (disperse_reload .AND. disperse_mindist >= 0.0) then
                    cand_info(ipos, iload)%disperse_dist = disperse_dist(iass)
                end if
                WRITE(NOUT, '(I5, 8ES13.5, 2I5, 3ES13.5)') iass, bocpwr(iass), eocpwr(iass), asspwr(iass,iload), pekpwr(iass, iload), pekass(iass, iload), &
                                                                bocden(iass), eocden(iass), pekden(iass,iload), pekpos(:,iass,iload), assbrn(iass), cand_info(ipos, iload)%rmsdif, &
                                                                pekfac(iass,iload)
            end do    
            num_candidate = 0
            margin_weight = margin_power_limit
            IF (objective_func .EQ. 'reactivity' .AND. .NOT. dif3d_check ) THEN
                write(NOUT,'(" margin for power limit =", ES13.5)') margin_weight
                npassed_power = 0
                npassed_burnup = 0
                npassed_distance = 0
                do ipos = 1, initial_pool_size
                    IASS = initial_pool(ipos)
                    tempcand => cand_info(iass2icand(IASS), iload)
                    lpass_all = .TRUE.
                    ! screen with peak power constraint
                    lpass_power = pass_fail(nout, tempcand%pekden, peak_power_limit*(1.0-margin_weight), eps_power, 1) .AND. &
                                pass_fail(nout, tempcand%pekpwr, assembly_power_limit*(1.0-margin_weight), eps_power, 1)
                    IF (.NOT. lpass_power) CYCLE
                    npassed_power = npassed_power + 1
                    passlist_power(npassed_power) = iass
                END DO
                IF (npassed_power .EQ. 0) THEN
                    objective_func = 'power_shape'
                    write(NOUT,'(" objective function change (reactivity -> power_shape) due to no candidate satisfies power margin")') 
                    write(NOUT,'(" and power limit is increased by 10% ")') 
                    peak_power_limit = peak_power_limit * 1.1
                    assembly_power_limit = assembly_power_limit * 1.1
                END IF
                npassed_power = 0
                npassed_burnup = 0
                npassed_distance = 0
                do ipos = 1, initial_pool_size
                    ! screen out positions inducing small power increase with refueling (indicating small reactivity addition)
                    ! check on this, 150.0 for 0.5 days, 210.0 for 1.0 day  -- tentative threshold, should be determined from history later
                    IASS = initial_pool(ipos)
                    tempcand => cand_info(iass2icand(IASS), iload)
                    lpass_all = .TRUE.
                    ! screen with peak power constraint
                    lpass_power = pass_fail(nout, tempcand%pekden, peak_power_limit*(1.0-margin_weight), eps_power, 1) .AND. &
                                pass_fail(nout, tempcand%pekpwr, assembly_power_limit*(1.0-margin_weight), eps_power, 1)
                    lpass_all = lpass_all .AND. lpass_power
                    IF (.NOT. lpass_power) CYCLE
                    npassed_power = npassed_power + 1
                    passlist_power(npassed_power) = iass
                    ! apply discharge burnup constraint
                    lpass_burnup = .TRUE.
                    if (consider_burnup .AND. .NOT. lrelease_burnup) then
                        lpass_burnup = pass_fail(nout, tempcand%assbrn, min_burnup, eps_power, 2)
                        lpass_all = lpass_all .AND. lpass_burnup
                    end if
                    IF (.NOT. lpass_all) CYCLE
                    npassed_burnup = npassed_burnup + 1
                    passlist_burnup(npassed_burnup) = iass
                    ! apply disperse distance
                    IF (disperse_reload) then
                        if(disperse_mindist < 0.0) then
                            !write(NOUT,'(" will exclude these from candidate refueling positions")')
                            lpass_distance = .NOT. lReaptedASM(iass, previous_pos, num_previous_site)
                        else!
                            lpass_distance = pass_fail(nout, tempcand%disperse_dist, disperse_mindist, eps_dist, 2)
                        end if
                        lpass_all = lpass_all .AND. lpass_distance
                    end if
                    IF (.NOT. lpass_all) CYCLE                
                    npassed_distance = npassed_distance + 1
                    passlist_distance(npassed_distance) = iass
                    ! store candidate that passed all the constraints except for keff
                    num_candidate = num_candidate + 1 
                    candidate_pos(num_candidate) = IASS
                    tempcand%lpassed = .TRUE.
                end do
                if (npassed_power>0) write(nout, '(a, I5, /, 1000I5)') ' Passed assembly (power)   ',npassed_power, passlist_power(1:npassed_power)
                if (npassed_burnup>0) write(nout, '(a, I5, /, 1000I5)') ' Passed assembly (burnup)  ',npassed_burnup, passlist_burnup(1:npassed_burnup)
                if (npassed_distance>0) write(nout, '(a, I5, /, 1000I5)') ' Passed assembly (distance)',npassed_distance, passlist_distance(1:npassed_distance)
            ELSE
                DO WHILE (num_candidate==0)
                    write(NOUT,'(" margin for power limit =", ES13.5)') margin_weight
                    npassed_power = 0
                    npassed_burnup = 0
                    npassed_distance = 0
                    do ipos = 1, initial_pool_size
                        IASS = initial_pool(ipos)
                        tempcand => cand_info(iass2icand(IASS), iload)
                        lpass_all = .TRUE.
                        ! screen with peak power constraint
                        lpass_power = pass_fail(nout, tempcand%pekden, peak_power_limit*(1.0-margin_weight), eps_power, 1) .AND. &
                                    pass_fail(nout, tempcand%pekpwr, assembly_power_limit*(1.0-margin_weight), eps_power, 1)
                        lpass_all = lpass_all .AND. lpass_power
                        IF (.NOT. lpass_all) CYCLE
                        npassed_power = npassed_power + 1
                        passlist_power(npassed_power) = iass
                        ! screen out positions inducing small power increase with refueling (indicating small reactivity addition)
                        ! check on this, 150.0 for 0.5 days, 210.0 for 1.0 day  -- tentative threshold, should be determined from history later
                        if (objective_func .NE. 'reactivity' .AND. objective_func .NE. 'reactivity_slow') THEN
                            lpass_minpower = (tempcand%asspwr - tempcand%eocpwr) >= 210.0
                            lpass_all = lpass_all .AND. lpass_minpower
                            IF (.NOT. lpass_all) CYCLE
                        END IF
                        ! apply discharge burnup constraint
                        lpass_burnup = .TRUE.
                        if (consider_burnup .AND. .NOT. lrelease_burnup) then
                            lpass_burnup = pass_fail(nout, tempcand%assbrn, min_burnup, eps_power, 2)
                            lpass_all = lpass_all .AND. lpass_burnup
                        end if
                        IF (.NOT. lpass_all) CYCLE
                        npassed_burnup = npassed_burnup + 1
                        passlist_burnup(npassed_burnup) = iass
                        ! apply disperse distance
                        IF (disperse_reload) then
                            if(disperse_mindist < 0.0) then
                                !write(NOUT,'(" will exclude these from candidate refueling positions")')
                                lpass_distance = .NOT. lReaptedASM(iass, previous_pos, num_previous_site)
                            else!
                                lpass_distance = pass_fail(nout, tempcand%disperse_dist, disperse_mindist, eps_dist, 2)
                            end if
                            lpass_all = lpass_all .AND. lpass_distance
                        end if
                        IF (.NOT. lpass_all) CYCLE                
                        npassed_distance = npassed_distance + 1
                        passlist_distance(npassed_distance) = iass
                        ! store candidate that passed all the constraints except for keff
                        num_candidate = num_candidate + 1 
                        candidate_pos(num_candidate) = IASS
                        tempcand%lpassed = .TRUE.
                    end do
                    if (npassed_power>0) write(nout, '(a, I5, /, 1000I5)') ' Passed assembly (power)   ',npassed_power, passlist_power(1:npassed_power)
                    if (npassed_burnup>0) write(nout, '(a, I5, /, 1000I5)') ' Passed assembly (burnup)  ',npassed_burnup, passlist_burnup(1:npassed_burnup)
                    if (npassed_distance>0) write(nout, '(a, I5, /, 1000I5)') ' Passed assembly (distance)',npassed_distance, passlist_distance(1:npassed_distance)
                    margin_weight = margin_weight * 0.5
                    IF (margin_weight<1E-10) EXIT
                    IF (margin_weight<1E-5) margin_weight = 0.0
                END DO
            END IF
            !=============================================== end of first exam
        end if ! manual or not
        if (num_candidate == 0) then
            search_failed = .TRUE.
            go to 900
        end if
        ! ========================== select candidate  (second exam) ============================================================
        !lselected = .FALSE. ! initialized at the beginning of Do loop for iload
        linit_persent = .TRUE.
        npersent_out_stack = 0
        ndif3d_out_stack = 0
        nsearch_iter = 0
        DO WHILE (.NOT. lselected)
            CALL SortCandidates(objective_func, peaks(1:num_candidate), candidate_pos, num_candidate, cand_info, iload)
            IF (linit_persent) write(NOUT,'(/," Potential candidates for PERSENT/DIF3D searches:",/,1000I5)'), candidate_pos(1:num_candidate)
            nsearch_iter= nsearch_iter+1
            WRITE(NOUT, '(/,A, I4, "...")') 'Iteration ', nsearch_iter
            if(consider_keff) then
                if(use_persent) then
                    ! work with PERSENT
                    ! buf_persent will be used to store the remaining lists of the potential candidates.
                    ! at each persent calculation, reactivity perturbation of 'num_persent_case' candidates are estimated and
                    ! remaining candidates are stored in candidates
                    IF (linit_persent) THEN
                        num_buf_persent = num_candidate
                        buf_persent(1:num_buf_persent)= candidate_pos(1:num_candidate)
                        linit_persent = .FALSE.
                    END IF
                    IF (num_buf_persent == 0 .AND. .NOT. lselected) THEN
                        search_failed = .true.
                        EXIT
                    END IF
                    done_with_persent = .false.
                    persent_time(1:3) = 0.0d0
                    do while (.not. done_with_persent)
                        if(num_persent_case > 0 .and. num_persent_case < num_buf_persent) then
                            ! shrink candidate pool and save untested cases
                            num_candidate = num_buf_persent - num_persent_case
                            candidate_pos(1:num_persent_case) = buf_persent(1:num_persent_case)
                            buf_persent(1:num_candidate) = buf_persent(num_persent_case+1:num_buf_persent)
                            buf_persent(num_candidate+1:num_buf_persent) = 0
                            num_buf_persent = num_candidate
                            num_candidate = num_persent_case
                        else
                            num_candidate = num_buf_persent
                            candidate_pos(1:num_candidate) = buf_persent(1:num_candidate)
                            num_buf_persent = 0
                            done_with_persent = .true.
                        end if
                        IF (dif3d_check) THEN
                            write(NOUT,'("  PERSENT/DIF3D check will be performed for the following candidates:",/,2x,1000I5)') candidate_pos(1:num_candidate)
                        ELSE
                            write(NOUT,'("  PERSENT check will be performed for the following candidates:",/,2x,1000I5)') candidate_pos(1:num_candidate)
                        END IF
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
                            npersent_out_stack = npersent_out_stack + 1
                            if(num_hex_sector > 0) then
                                write(persent_out_stack(npersent_out_stack),'(I6,3X,2I6,4X,F13.5,5X,F13.5,4X,F13.5)') IASS, hexagon_rp(:,IASS),fdum,fdum2,deltak(IASS,iload)
                            else
                                call reload_ass2xy(IASS, num_xnode, IX, IY)
                                write(persent_out_stack(npersent_out_stack),'(I6,3X,2I6,4X,F13.5,5X,F13.5,4X,F13.5)') IASS, IX, IY, fdum,fdum2,deltak(IASS,iload)
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
                        else if(num_candidate == 0 .and. num_buf_persent == 0) then
                            search_failed = .true.
                            done_with_persent = .true.
                            EXIT
                        end if
                    end do
                    call system_command(NOUT,'mv  user.NAFLUX  init.NAFLUX')   ! saved as initial guess for next cycle 
                    call system_command(NOUT,'rm  user.NHFLUX')   ! will have new NHFLUX in next cycle 
                    IF (search_failed) THEN
                        write(NOUT, '(/,A, I6)') 'All the candidates cannot meet reactivity constraints from the PERSENT calculations!'
                        EXIT
                    END IF
                else ! if persent is not used...
                    !if persent is not used, all the candidates passed the first exam
                    write(NOUT, '(/,A, I6)') 'Caution !!! all the candidates will be screened by dif3d calculation'
                end if
            end if
            
            write(nout, *)
            ! ====================== determine final position of refueling based on objective_func
            if(.not. dif3d_check) then
                WRITE(NOUT,'(/,A,1000i6)') 'Select refueling assembly among (w/o DIF3D check) :', candidate_pos(1:num_candidate)
                reload_pos = select_candidate(nout, candidate_pos, num_candidate, cand_info, iload)
                peak_pwrden = pekden(reload_pos,iload)
                peak_asspwr = pekpwr(reload_pos,iload)
                reload_pekass = pekass(reload_pos,iload)
                reload_pekpos = pekpos(1:2,reload_pos,iload)
                if(consider_keff .and. use_persent) peak_keff = eoc_keff + deltak(reload_pos,iload)
                lselected = .true.
            else         ! =================== final check with dif3d calculation
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
                if(objective_func == 'power_shape') then
                    write(NOUT,100) 'Warning: DIF3D check for objective function <power_shape> not completed, &
                    &relative power change from reference based on DIF3D results not added'
                    call abort 
                end if
                CALL SortCandidates(objective_func, peaks(1:num_candidate), candidate_pos, num_candidate, cand_info, iload)

                WRITE(NOUT,'(/,A, 100I6)') 'Sort candidate_pos : ', candidate_pos(1:num_candidate)
                ! run dif3d with replaced assembly at candidate positions (from most to least probable), stop when all constraints are ensured
                ! now use scratch_pool to record peak assembly position for each reloading
                call reload_rundif3d(NOUT,num_candidate,candidate_pos,iload,RFUEL(iload),keffs(1:num_candidate),peaks(1:num_candidate),&
                &   scratch_pool(1:num_candidate), peakden(1:num_candidate), scratch_pool2(:,1:num_candidate), num_scratch)    
                !call system_command(NOUT, 'rm -r '//dif3d_dir)      ! remove this if like to see DIF3D output

                WRITE(NOUT,'(/,A, 100I6)') 'reload candidate_pos : ', candidate_pos(1:num_candidate)
                num_candidate = num_scratch
                peak_asspwr = peaks(1) * 1.5
                peak_pwrden = peakden(1) * 1.5    ! make sure not miss the first one
                peak_burnup = -1.0
                peak_keff   = keffs(1) * 0.5      ! make sure not miss the first one
                peak_pekfac = pekfac(abs(candidate_pos(ipos)),iload)
                do ipos = 1, num_candidate
                    IASS = candidate_pos(ipos)
                    if(IASS > 0) then
                        cdum = 'Yes'
                    else
                        cdum = 'No'
                        IASS = -IASS
                    end if
                    ndif3d_out_stack = ndif3d_out_stack + 1
                    if(num_hex_sector > 0) then
                        write(dif3d_out_stack(ndif3d_out_stack),'(I5,2I6,3X,F11.5,3X,ES13.5,3X,I6,13X,ES13.5,2I10,5X,A)') IASS,hexagon_rp(1:2,IASS),keffs(ipos)-eoc_keff,&
                        &   peaks(ipos),scratch_pool(ipos), peakden(ipos), scratch_pool2(:,ipos), cdum
                    else
                        call reload_ass2xy(IASS, num_xnode, IX, IY)
                        write(dif3d_out_stack(ndif3d_out_stack),'(I5,2I6,3X,F11.5,3X,ES13.5,3X,I6,13X,ES13.5,2I10,5X,A)') IASS, IX, IY, keffs(ipos)-eoc_keff, &
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

                        ! to maximize reactivity
                        if(objective_func == 'reactivity_slow') then
                            fdum = (keffs(ipos) - peak_keff) - conv_peak_rho*wgt_pekfac*(pekfac(IASS, iload)-peak_pekfac)
                            if(fdum > eps_keff) then
                                update_position = .true.
                            else if(abs(fdum) <= eps_keff) then  ! peak candidate with smaller peak
                                fdum = peakden(ipos) / peak_pwrden - wgt_pekfac*(pekfac(IASS, iload)-peak_pekfac)
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
                if(num_scratch == 0) THEN
                    IF (num_buf_persent == 0) then
                        search_failed = .true.
                        WRITE(NOUT,'(A)') '     PERSENT/DIF3D Search FAIL and there is no candidate...'
                        EXIT
                    ELSE
                        WRITE(NOUT,'(A)') '     PERSENT/DIF3D Search FAIL...'
                    END IF
                else
                    lselected = .true.
                    WRITE(NOUT,'(A)') ' PERSENT/DIF3D Search SUCCESS!'
                end if
            end if
        end DO

        IF (npersent_out_stack>0) THEN
            ! print persent table header
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
            DO ipos = 1, npersent_out_stack
                write(nout, '(A)') TRIM(persent_out_stack(ipos))
            END DO
        END IF

        IF (ndif3d_out_stack>0) THEN
            ! print dif3d table header
            write(NOUT,'(/,20X,"*** DIF3D calculations with reloaded assembly at candidate positions ***")')
            if(num_hex_sector > 0) then
                write(NOUT,'(/,36X,"Peak assembly integrated power     Peak node-averaged power density    Valid",/,&
            &"Index  IRING  IPOS    Delta_Keff    Value (W)   Assembly Location      Value (W/cc)    Assembly    Node    Candidate?",/)')
            else
                write(NOUT,'(/,36X,"Peak assembly integrated power     Peak node-averaged power density    Valid",/,&
            &"Index    IX    IY     Delta_Keff    Value (W)   Assembly Location      Value (W/cc)    Assembly    Node    Candidate?",/)')
            end if
            DO ipos = 1, ndif3d_out_stack
                write(nout, '(A)') TRIM(dif3d_out_stack(ipos))
            END DO
        END IF
        
900     continue

        ! IF (search_failed .AND. objective_func == 'reactivity') THEN

        !     search_failed = .TRUE.
        ! END IF


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
            else if(objective_func == 'reactivity_slow') then
                write(NOUT,'(" reloading assembly position is chosen to have maximum reactivity increase with minimum peaking factor increase")')
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
    end do ! iload 

    ! place holder for now (pick one among multiple reload composition and position)
    write(NOUT,'(/,20X,"*** Summary of reloaded assembly ***",/)')
    write(NOUT,'("reloaded         reload_position     reload_regions",/,"composition      IASS   IX   IY      names")')
    do iload = 1, num_reload_comp
        if(iSearch(iload) > 0) &
        &   write(NOUT,'(I2,3X,A8,3X,3I5,6X,100(A8,1X))') iload, RFUEL(iload)%name, reload_pos, IX, IY, reload_reg(1:nreg)
    end do
    
  CONTAINS
        function Select_Candidate(nout, candidate_pos, num_candidate, cand_info, iload)
            integer :: nout
            integer :: num_candidate, iload
            integer :: candidate_pos(num_candidate)
            TYPE(candidate_info_type), pointer :: cand_info(:,:)
            integer :: Select_Candidate
            ! local
            integer :: ipos, iass, icand
            real(8)  :: fdum, fdum2
            !
            Select_Candidate = candidate_pos(1)
            icand = iass2icand(Select_Candidate)
            select case (trim(objective_func))
            case ('assembly_power')
                peak_asspwr= cand_info(icand, iload)%pekpwr
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
                    icand = iass2icand(iass)
                    fdum = cand_info(icand, iload)%pekpwr / peak_asspwr - 1.0D0
                    if(fdum < -eps_power) then
                        Select_Candidate = IASS
                        peak_asspwr = pekpwr(IASS,iload)
                    end if
                end do
            case ('peak_power')
                peak_pwrden = pekden(Select_Candidate,iload)
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
                    icand = iass2icand(iass)
                    fdum = pekden(IASS,iload) / peak_pwrden - 1.0D0
                    if(fdum < -eps_power) then
                        Select_Candidate = IASS
                        peak_pwrden = pekden(IASS,iload)
                    end if
                end do
            case ('burnup')
                peak_burnup = assbrn(Select_Candidate)
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
          icand = iass2icand(iass)
                    fdum = assbrn(IASS) / peak_burnup - 1.0D0
                    if(fdum > eps_power) then
                        Select_Candidate = IASS
                        peak_burnup = assbrn(IASS)
                    end if
                end do
            case ('reactivity')
                peak_keff = deltak(Select_Candidate,iload)
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
                    icand = iass2icand(iass)
                    if(deltak(IASS,iload) - peak_keff > eps_keff) then
                        Select_Candidate = IASS
                        peak_keff = deltak(IASS,iload)
                    else if(abs(deltak(IASS,iload) - peak_keff) <= eps_keff) then
                        if(assbrn(IASS)/assbrn(Select_Candidate) - 1.0d0 > eps_power) then
                            Select_Candidate = IASS
                            peak_keff = deltak(IASS,iload)
                        end if
                    end if
                end do
            case ('reactivity_slow')
                peak_keff = deltak(Select_Candidate,iload)
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
                    icand = iass2icand(iass)
                    fdum = (deltak(IASS,iload) - peak_keff) - (pekfac(iass,iload)-pekfac(Select_Candidate,iload))*conv_peak_rho
                    fdum2 = ABS(deltak(IASS,iload) - peak_keff) - (pekfac(iass,iload)-pekfac(Select_Candidate,iload))*conv_peak_rho
                    if(fdum > eps_keff) then
                        Select_Candidate = IASS
                        peak_keff = deltak(IASS,iload)
                    else if(fdum2 <= eps_keff) then
                        if(assbrn(IASS)/assbrn(Select_Candidate) - 1.0d0 > eps_power) then
                            Select_Candidate = IASS
                            peak_keff = deltak(IASS,iload)
                        end if
                    end if
                end do
            case ('reactivity_rev')
                peak_keff = deltak(Select_Candidate,iload)
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
                    icand = iass2icand(iass)
                    fdum = (deltak(IASS,iload) - peak_keff) - (pekfac(iass,iload)-pekfac(Select_Candidate,iload))*conv_peak_rho
                    fdum2 = ABS(deltak(IASS,iload) - peak_keff) - (pekfac(iass,iload)-pekfac(Select_Candidate,iload))*conv_peak_rho
                    if(fdum > eps_keff) then
                        Select_Candidate = IASS
                        peak_keff = deltak(IASS,iload)
                    else if(fdum2 <= eps_keff) then
                        if(assbrn(IASS)/assbrn(Select_Candidate) - 1.0d0 > eps_power) then
                            Select_Candidate = IASS
                            peak_keff = deltak(IASS,iload)
                        end if
                    end if
                end do
            case ('power_shape')
                fdum = cand_info(icand, iload)%rmsdif
                fdum2 = cand_info(icand, iload)%pekpwr
                update_position = .false.
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
          icand = iass2icand(iass)
                    if((fdum - cand_info(icand, iload)%rmsdif) > eps_power_diff) then 
                        update_position = .true.
                    else if(abs(fdum - cand_info(icand, iload)%rmsdif) < eps_power_diff) then
                        if(fdum2 / pekpwr(IASS,iload) - 1.0d0 > eps_power) update_position = .true.
                    end if
                    if(update_position) then
                        Select_Candidate = IASS
                        fdum = cand_info(icand, iload)%rmsdif
                        fdum2 = pekpwr(IASS,iload)
                    end if
                end do
            case ('power_deep')     ! Moltex suggested option
                fdum2 = eocpwr(Select_Candidate) / refpwr(Select_Candidate)
                do ipos = 2, num_candidate
                    IASS = candidate_pos(ipos)
          icand = iass2icand(iass)
                    fdum = eocpwr(IASS) / refpwr(IASS)
                    if(fdum - fdum2 < -eps_power_diff) then
                        Select_Candidate = IASS
                        fdum2 = fdum
                    end if
                end do
            case default
                write(NOUT,'(A)') 'Error: unrecognized objective_func '//objective_func
                call abort('')
      end select
        END function
    
        SUBROUTINE SortCandidates(objective_func, peaks, candidate_pos, num_candidate, cand_info, iload)
            character(len=*) :: objective_func  ! 'ascend' / 'descend'
            integer :: num_candidate, iload
            real(8) :: peaks(num_candidate)
            integer :: candidate_pos(num_candidate)
            TYPE(candidate_info_TYPE), POINTER :: cand_info(:, :)
            ! local 
            integer :: ipos, iass, icand
            if (num_candidate==1) return
            select case (trim(objective_func))
                case ('assembly_power') ! sort candidate positions in order of increasing peak assembly power after refueling for batch wise run
                    do ipos = 1, num_candidate
            icand = iass2icand(candidate_pos(ipos))
                        peaks(ipos) = cand_info(icand, iload)%pekpwr
                    end do
                    call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend', eps_power)
                case ('peak_power') ! sort candidate positions in order of increasing peak power density after refueling for batch wise run
                    do ipos = 1, num_candidate
            icand = iass2icand(candidate_pos(ipos))
                        peaks(ipos) = cand_info(icand, iload)%pekden
                    end do
                    call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend', eps_power)
                case ('burnup') ! sort candidate positions in order of decreasing discharge burnup for batch wise run
                    do ipos = 1, num_candidate
            icand = iass2icand(candidate_pos(ipos))
                        peaks(ipos) = cand_info(icand, iload)%assbrn
                    end do
                    call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'descend', eps_power)
                case ('reactivity') ! sort candidate positions in order of decreasing relative change of assembly powers due to refueling for batch wise run
                    do ipos = 1, num_candidate
            icand = iass2icand(candidate_pos(ipos))
                        peaks(ipos) = cand_info(icand, iload)%asspwr / cand_info(icand, iload)%eocpwr - 1.0d0     ! relative power change due to refueling, indicator of added reactivity
                    end do
                    call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'descend', -eps_keff)
                case ('reactivity_slow') ! sort candidate positions in order of decreasing relative change of assembly powers due to refueling for batch wise run
                    do ipos = 1, num_candidate
            icand = iass2icand(candidate_pos(ipos))
                        peaks(ipos) = cand_info(icand, iload)%asspwr / cand_info(icand, iload)%eocpwr - wgt_pekfac * (cand_info(icand, iload)%pekfac - eocpekfac)
                    end do
                    call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'descend', -eps_keff)
                case ('power_shape') ! sort candidate positions in order of increasing deviation of power distribution
                    avgpwr = PWDINT%POWER / num_potential_site  ! approximate average fuel assembly power
                    digit_shift = 10**(-floor(log10(eps_power_diff)))      ! reasonable assumption: eps_power_diff < 1.0
                    do ipos = 1, num_candidate
            icand = iass2icand(candidate_pos(ipos))
                        IASS = candidate_pos(ipos)
                        call reload_power_change(iload, IASS, cand_info(icand, iload)%maxup, cand_info(icand, iload)%maxpos, cand_info(icand, iload)%rmsdif)
                        peaks(ipos) = floor(cand_info(icand, iload)%rmsdif*digit_shift+0.5d0)*10 + (cand_info(icand, iload)%pekpwr / avgpwr) ! the latter should be well below 10.0
                    end do
                    call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend', -eps_power)
                case ('power_deep')  ! sort candidate positions in order of increasing relative assembly power (Power_eoc/Power_ref) for batch wise run
                    do ipos = 1, num_candidate
            icand = iass2icand(candidate_pos(ipos))
                        peaks(ipos) = cand_info(icand, iload)%eocpwr / cand_info(icand, iload)%refpwr - 1.0d0 ! power change relative to reference power
                    end do
                    call insertSort(num_candidate, peaks(1:num_candidate), candidate_pos(1:num_candidate), 'ascend', -eps_power_diff)
                case default
                    write(NOUT,'(a)') 'Note: batch wise PERSENT run not implemented for this objective_func,&
                    & will estimate reactivity changes in ONE PERSENT run'
                    num_persent_case = num_candidate
                end select
        END SUBROUTINE
        FUNCTION pass_fail(nout, var, constraint, tol, mode)
            INTEGER :: nout        ! device number for output
            real(8) :: var         ! variable to examine pass or fail
            real(8) :: constraint  ! constraint of variable
            real(8) :: tol         ! tolerance of variable
            integer :: mode        ! 1: pass if var < constraint / 2: pass if var > constraint 
            LOGICAL :: pass_fail
            
            real(8) :: fdum
            !
            pass_fail = .FALSE.
            if (mode == 1) then
                fdum = var / constraint - 1.0D0
                if(fdum <= tol) pass_fail = .TRUE.
            else if (mode == 2) then
                fdum = var / constraint - 1.0D0
                if(fdum >= -tol) pass_fail = .TRUE.
            else
                write(100, '(a, i)') 'pass_fail - unsupported mode : ', mode
                call abort          
            end if
            RETURN
        End FUNCTION
        !
        FUNCTION lReaptedASM(IASS, previous_pos, num_previous_site)
            INTEGER :: IASS, num_previous_site
            INTEGER, POINTER :: previous_pos(:)
            LOGICAL :: lReaptedASM
            INTEGER :: kk
            lReaptedASM = .FALSE.
            do kk = 1, num_previous_site
                if(IASS == previous_pos(kk)) then
                    lReaptedASM = .TRUE.
                    RETURN
                end if
            end do
        END FUNCTION
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

End Subroutine

end module