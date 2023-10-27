Subroutine reload_potential_site(NOUT)
    ! locate potential reloading site, i.e., fuel assemblies, using information from input and GEODST
    ! make sure fuel assemblies are modeled individually, no region is assigned to multiple assemblies
    use m_reload
    implicit none
    integer :: NOUT
    integer :: IZ, IX, IY, IREG, ii, jj
    character(len=8) :: ANAME
    logical :: yes
    integer, allocatable :: tmp(:)

    100 format('[reload_potential_site]...',A,1X,2I6)
    num_xnode = GEODST%NINTI
    num_ynode = GEODST%NINTJ
    num_xynod = num_xnode * num_ynode
    allocate(assembly_map(num_xnode, num_ynode), tmp(num_xynod))

    assembly_map(:,:) = 0
    bot_active_node = 0
    do IZ = 1, GEODST%NINTK
        do IY = 1, num_ynode
            do IX = 1, num_xnode
                IREG = GEODST%MR(IX,IY,IZ)
                if(IREG == 0) cycle     ! out of core domain
                if(IZ == 1) assembly_map(IX,IY) = -1    ! mark other assemblies 
                ANAME = LABELS%REGNAM(IREG)
                do ii = 1, num_fuel_region
                    if(ANAME == fuel_regions(ii)) exit
                end do
                if (ii <= num_fuel_region) then
                    assembly_map(IX,IY) = reload_type   ! mark fuel assembly
                    if(bot_active_node == 0) bot_active_node = IZ
                    top_active_node = IZ
                end if
            end do
        end do
    end do

    ! remove symmetric sites in hexagonal geometry and record all fuel assembly positions
    ii = 0
    num_potential_site = 0
    do IY = 1, num_ynode 
        do IX = 1, num_xnode 
            ii = ii + 1
            if(hexagon_rp(1, ii)==0) assembly_map(IX,IY) = 0
            if(assembly_map(IX,IY) == reload_type) then
                num_potential_site = num_potential_site + 1
                tmp(num_potential_site) = ii
            end if
        end do
    end do
    allocate(fuel_pos(num_potential_site))
    fuel_pos = tmp(1:num_potential_site)
    call reload_candidate_reorder(num_potential_site,fuel_pos)

    ! check whether individual fuel assemblies are modeled
    deallocate(tmp)
    allocate(tmp(GEODST%NREG))
    tmp(:) = 0
    jj = 0
    do ii = 1, num_potential_site
        call reload_ass2xy(fuel_pos(ii), GEODST%NINTI, IX, IY)
        do IZ = bot_active_node, top_active_node
            IREG = GEODST%MR(IX,IY,IZ)
            if(IZ > bot_active_node) THEN
              IF (IREG == GEODST%MR(IX,IY,IZ-1)) cycle
            END IF
            tmp(IREG) = tmp(IREG) + 1
            if(tmp(IREG) > 1) then
                write(NOUT,100) 'Error: fuel region '//LABELS%REGNAM(IREG)//' is assigned to multiple assemblies'
                jj = 1
            end if
        end do
    end do
    if(jj > 0) then
        write(NOUT,100) 'Note: to compute assembly-wise discharge burnups and to replace single assembly at one time, &
        &a fuel region cannot be assigned to multiple fuel assemblies'
        call abort
    end if

    ! miscellaneous checks
    yes = .false.
    do ii = 1, num_previous_site
        jj = previous_pos(ii)
        call reload_ass2xy(jj, num_xnode, IX, IY)
        if(assembly_map(IX,IY) /= reload_type) then
            if(num_hex_sector == 0) then
                write(NOUT,100) 'Error: fresh fuel was reloaded at non-fuel assembly position, (IX, IY) =', IX, IY
            else
                write(NOUT,100) 'Error: fresh fuel was reloaded at non-fuel assembly position, (IRING, IPOS) =', hexagon_rp(1:2,jj)
            end if
            yes = .true.
        end if
    end do
    if(manual_reload) then
        if(num_hex_sector > 0) then
            call GetIJPOS(manual_reload_pos(1),manual_reload_pos(2),IX,IY,num_hex_sector,num_xnode,num_ynode,NOUT)
            if(IX < 0 .or. IY < 0) then
                write(NOUT,100) 'Error: user specified refueling position is out-of-domain, (IRING, IPOS) =', manual_reload_pos
                yes = .true.
            else if(assembly_map(IX,IY) /= reload_type) then
                write(NOUT,100) 'Error: fresh fuel is manually reloaded at non-fuel assembly position, (IRING, IPOS) =',manual_reload_pos
                yes = .true.
            end if
            manual_reload_pos(1) = IX
            manual_reload_pos(2) = IY
        else
            IX = manual_reload_pos(1)
            IY = manual_reload_pos(2)
            if(assembly_map(IX,IY) /= reload_type) then
                write(NOUT,100) 'Error: fresh fuel is manually set to be reloaded at non-fuel assembly position, (IX, IY) =', IX, IY
                yes = .true.
            end if
        end if
    end if

    ! check prescribed candidate positions and convert (IRING, IPOS) to (IX,IY)
    if(candidate_range_size > 0) then   ! user specifies candidate range
        if(num_hex_sector > 0) then
            do ii = 1, candidate_range_size
                call GetIJPOS(candidate_range(1,ii),candidate_range(2,ii),IX,IY,num_hex_sector,num_xnode,num_ynode,NOUT)
                if(IX < 0 .or. IY < 0) then
                    write(NOUT,100) 'Error: user specified candidate position is out-of-domain, (IRING, IPOS) =', candidate_range(:,ii)
                    yes = .true.
                else if(assembly_map(IX,IY) /= reload_type) then
                    write(NOUT,100) 'Error: user specified candidate positions include a non-fuel assembly, (IRING, IPOS) =', candidate_range(:,ii)
                    yes = .true.
                end if
                candidate_range(1,ii) = IX
                candidate_range(2,ii) = IY
            end do
        else
            do ii = 1, candidate_range_size
                IX = candidate_range(1,ii)
                IY = candidate_range(2,ii)
                if(assembly_map(IX,IY) /= reload_type) then
                    write(NOUT,100) 'Error: user specified candidate positions include a non-fuel assembly, (IX, IY) =', IX, IY
                    yes = .true.
                end if
            end do
        end if

        allocate(candidate_range_pos(candidate_range_size))
        do ii = 1, candidate_range_size
            candidate_range_pos(ii) = (candidate_range(2,ii)-1) * num_xnode + candidate_range(1,ii)
        end do
    end if

    if(yes) call abort

End Subroutine

! ------------------------------------------------------------------------------
! reorder candidate positions in the order of increasing ring number and 
! increasing position number in each ring
Subroutine reload_candidate_reorder(num_candidate, candidate_pool)
    use m_utility
    use m_reload
    implicit none
    integer :: num_candidate, candidate_pool(num_candidate)
    ! local
    integer :: IR, IP, IASS, ii, digit_shift
    real(8) :: sorting(num_candidate)

    IP = maxval(hexagon_rp(2,:))
    digit_shift = 1000
    if(IP > 1000) digit_shift = 10000

    do ii = 1, num_candidate
        IASS = candidate_pool(ii)
        IR = hexagon_rp(1,IASS)
        IP = hexagon_rp(2,IASS)
        sorting(ii) = 1.0*(IR*digit_shift + IP)
    end do
    call insertSort(num_candidate, sorting, candidate_pool, 'ascend', -0.5d0)

End Subroutine