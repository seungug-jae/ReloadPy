Program reload
    use m_reload
    use GEODST_IO
    use RTFLUX_IO
    use NHFLUX_IO
    use LABELS_IO
    use PWDINT_IO
    use find_position
    implicit none 
    integer :: NOUT
    data NOUT /11/
    logical :: yes
    real(8) :: tic, toc

    100 format('[reload]...',A,1X,I6)

    call CPU_TIME(tic)
    open(NOUT,file='reload.out',status='replace',form='formatted')
    call reload_input(NOUT)
    if(null_refueling) then
        write(NOUT,'(" Constraints for reloading position: null_refueling")')
    else
        write(NOUT,'(" Constraints for reloading position: peak power limit, assembly power limit")',advance='NO')
        if(consider_burnup) write(NOUT,'(", discharge burnup")',advance='NO')
        if(consider_keff .or. consider_worth) write(NOUT,'(", reactivity")',advance='NO')
        if(adjust_cycle_length) write(NOUT,'(", cycle length")',advance='NO')
        if(disperse_reload) write(NOUT,'(", assembly distance")',advance='NO')
        if(manual_reload) write(NOUT,'(", user selection")', advance='NO')
    end if
    write(NOUT,'(/,"!!! Note that peak power density was determined as the maximum node-averaged quantity !!!")')
    write(NOUT,'("number of fresh fuel compositions =",I4,/,100(A,1X))') num_reload_comp, RFUEL(1:num_reload_comp)%name
    write(NOUT,'("primary objective: ",A)') objective_func
    write(NOUT,*)

    ! read GEODST
    call GEODST_IMPORT_BINARY(NOUT, GEODST, 'GEODST')
    if(GEODST%NINTI /= GEODST%NCINTI .or. GEODST%NINTJ /= GEODST%NCINTJ) then
        write(NOUT,100) 'Error: flux calculation is not performed with DIF3D-Nodal nor VARIANT option, mesh/assembly mapping is not implemented' 
        call abort
    end if
    call GEODST_NRASS1(GEODST)      ! convert region assignment map to node-wise 
    !call GEODST_ASSIGNPRINTINFO(NOUT)
    !call GEODST_PRINT(GEODST)

    ! compute macxs of reloaded composition
    ! may need compute COMPXS for reloaded fuel if perturbation method is used for evaluating reactivity change
    !if(consider_keff .and. (.not.use_persent)) then
    !    call reload_macxs(NOUT, 2)
    !else
        call reload_macxs(NOUT, 1)
    !end if

    ! read flux data
    if(flux_data == 'RTFLUX') then
        call RTFLUX_IMPORT_BINARY(NOUT, RTFLUX, 'RTFLUX')
        if(GEODST%NINTI /= RTFLUX%NINTI .or. GEODST%NINTJ /= RTFLUX%NINTJ .or. GEODST%NINTK /= RTFLUX%NINTK) then
            write(NOUT,100) 'Error: inconsistent geometry between GEODST and RTFLUX'
            write(NOUT,'("GEODST NINTI, NINTJ, NINTK =",3I4)') GEODST%NINTI, GEODST%NINTJ, GEODST%NINTK
            write(NOUT,'("RTFLUX NINTI, NINTJ, NINTK =",3I4)') RTFLUX%NINTI, RTFLUX%NINTJ, RTFLUX%NINTK
            call abort
        end if
        if(kerma%ngroup /= RTFLUX%NGROUP) then
            write(NOUT,100) 'Error: inconsistent group between ISOTXS and RTFLUX'
            write(NOUT,'("NGROUP in ISOTXS and RTFLUX are ",I4," and ",I4," respectively")') kerma%ngroup, RTFLUX%NGROUP
            call abort
        end if
        eoc_keff = RTFLUX%EFFK
        call RTFLUX_IMPORT_BINARY(NOUT, boc_RTFLUX, 'boc_RTFLUX')
        boc_keff = boc_RTFLUX%EFFK

    else if(flux_data == 'NHFLUX') then
        inquire(file='NHFLX0', exist=yes)
        if(yes) then
            call NHFLUX_IMPORT(NOUT, NHFLUX, 'NHFLX0', .true.)
        else
            call NHFLUX_IMPORT(NOUT, NHFLUX, 'NHFLUX', .true.)
        end if
        if(GEODST%NINTI /= NHFLUX%NINTI .or. GEODST%NINTJ /= NHFLUX%NINTJ .or. GEODST%NINTK /= NHFLUX%NINTK) then
            write(NOUT,100) 'Error: inconsistent geometry between GEODST and NHFLUX'
            call abort 
        end if
        if(kerma%ngroup /= NHFLUX%NGROUP) then
            write(NOUT,100) 'Error: inconsistent group between ISOTXS and NHFLUX'
            call abort
        end if
        eoc_keff = NHFLUX%EFFK
    end if

    ! read PWDINT
    call PWDINT_IMPORT_BINARY(NOUT,PWDINT,'PWDINT')
    if(GEODST%NINTI /= PWDINT%NINTI .or. GEODST%NINTJ /= PWDINT%NINTJ .or. GEODST%NINTK /= PWDINT%NINTK) then
        write(NOUT,100) 'Error: inconsistent geometry between GEODST and PWDINT'
        call abort 
    end if
    if(flux_data == 'RTFLUX' .and. RTFLUX%POWER /= PWDINT%POWER) then
        write(NOUT,100) 'Error: inconsistent nominal power between RTFLUX and PWDINT'
        call abort 
    end if
    if(flux_data == 'NHFLUX' .and. NHFLUX%POWER /= PWDINT%POWER) then
        write(NOUT,100) 'Error: inconsistent nominal power between RTFLUX and PWDINT'
        call abort 
    end if
    call PWDINT_IMPORT_BINARY(NOUT,boc_PWDINT,'boc_PWDINT')

    ! read LABELS for region names
    call LABELS_IMPORT(NOUT, LABELS, 'LABELS')

    ! make map for hexagonal assembly
    if(GEODST%IGOM == 10 .or. GEODST%IGOM == 18) then
        if(num_hex_sector == 0) then
            write(NOUT,100) 'Error: in hexagonal geometry, number of 60 degree sectors (num_hex_sector) have to be specified in &control'
            call abort
        end if
        allocate(hexagon_rp(2, GEODST%NCINTI*GEODST%NCINTJ))
        call MapIJ2RP(hexagon_rp, num_hex_sector, GEODST%NCINTI, GEODST%NCINTJ, NOUT)
    end if

    ! locate fuel assemblies
    call reload_potential_site(NOUT)

    ! compute total power of reloaded assembly
    call reload_power(NOUT)

    ! compute assembly discharge burnup
    call reload_burnup(NOUT)

    ! compute peak fast fluence (node level)
    call reload_fast_fluence(NOUT, RTFLUX%NINTI, RTFLUX%NINTJ, RTFLUX%NINTK, RTFLUX%NGROUP, boc_RTFLUX%FREG, RTFLUX%FREG)

    if(null_refueling) then
        call reload_null(NOUT)
    else ! determine reloading position 
        if(objective_func == 'reactivity') then
            call reload_position_2(NOUT)
        elseif (objective_func == 'reactivity_slow') then
            call reload_position_2(NOUT)
        else
            call reload_position(NOUT)
        end if
    end if

    call CPU_TIME(toc)
    write(NOUT,'(/,"total execution time (sec) =",ES13.5)') toc-tic

End Program