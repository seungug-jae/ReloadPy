Subroutine reload_persent_output(NOUT, filename, user_keff)
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
    
    100 format('[reload_persent_output]...',A,1X,I6)

    ! read PERSENT output for delta_rho and accumulate reactivity worth by replacing each fuel assembly
    UNIT = NE_Kind_GetFreeLogicalUnit()
    open(UNIT,file=filename, status='old',form='formatted')
    if(.not. allocated(deltak)) then
        ! if (measure_uncertainty) then
        !     if(consider_worth) then
        !         allocate(deltak(num_xynod, num_reload_comp+num_test_comp+1))
        !     else
        !         allocate(deltak(num_xynod, num_reload_comp+num_test_comp))
        !     end if
        ! else
            if(consider_worth) then
                allocate(deltak(num_xynod, num_reload_comp+1))
            else
                allocate(deltak(num_xynod, num_reload_comp))
            end if
        ! end if
        deltak(:,:) = 0.0d0
    end if

    keff = user_keff    ! if keff < 0.0,  will get it from PERSENT output
    do jj = 1, num_pert_cases
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
        deltak(IASS,iload) = keff * keff * drho / (1.0d0 - keff*drho)
    end do
    ! if (measure_uncertainty) then
    !     do jj = num_pert_cases+1, num_pert_cases + num_test_comp
    !         call SkipTo('[PERSENT]...Problem', .true., UNIT, IOS)
    !         call SkipTo('[PERSENT]...Performing the DIF3D-VARIANT numerator/denominator operations',.true., UNIT,IOS)
    !         read(UNIT,'(A)') line
    !         !ii = index(line,'is')
    !         !read(line(ii+2:line_len),*) drho
    !         ! seems fixed format is used in PERSENT output
    !         ii = index(line,'RF_FA_')
    !         read(line(ii+6:line_len),'(I5,A1,I2)') IASS, splitter, iload
    !         read(line(58:line_len),*) drho
    !         if(keff < 0.0) read(line(120:line_len),*) keff
    !         deltak(IASS,iload) = keff * keff * drho / (1.0d0 - keff*drho)
    !     end do
    ! end if
    close(UNIT)

End Subroutine


! ----------------------------------------------------------------------
Subroutine reload_persent_input(filename, iflux, num_pool, pool)
    ! make persent input for zone perturbation problem
    use m_reload
    use NE_Kind, only: NE_Kind_GetFreeLogicalUnit
    implicit none
    integer :: iflux      ! iflux = 0 / 1 / 2  prepare input with no flux / forward flux / both forward and adjoint flux provided
    character(len=*) :: filename
    integer :: num_pool, pool(num_pool)
    ! local
    integer :: UNIT, NIN, IDUM, IOS, NCMP, NCMP2
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

    NCMP = num_reload_comp
    if(consider_worth) NCMP = num_reload_comp + 1
    do iload = 1, NCMP
        ANAME = RFUEL(iload)%name
        jj = (RFUEL(iload)%niso + 9) / 10
        kk = 0
        do ii = 1, jj-1
            write(UNIT,'("NEW_ZONE  ",A,1X,10(1X,A8,ES13.5))') ANAME, ((RFUEL(iload)%isonam(kk+ij), RFUEL(iload)%density(kk+ij)),ij=1,10)
            kk = kk + 10
        end do
        if(kk < RFUEL(iload)%niso) &
        write(UNIT,'("NEW_ZONE  ",A,1X,10(1X,A8,ES13.5))') &
        &   ANAME, ((RFUEL(iload)%isonam(ij), RFUEL(iload)%density(ij)),ij=kk+1,RFUEL(iload)%niso)
        write(UNIT,*)
    end do
    ! if (measure_uncertainty) then
    !     NCMP2 = num_test_comp
    !     do iload = 1, NCMP2
    !         ANAME = TFUEL(iload)%name
    !         jj = (TFUEL(iload)%niso + 9) / 10
    !         kk = 0
    !         do ii = 1, jj-1
    !             write(UNIT,'("NEW_ZONE  ",A,1X,10(1X,A8,ES13.5))') ANAME, ((TFUEL(iload)%isonam(kk+ij), TFUEL(iload)%density(kk+ij)),ij=1,10)
    !             kk = kk + 10
    !         end do
    !         if(kk < TFUEL(iload)%niso) &
    !         write(UNIT,'("NEW_ZONE  ",A,1X,10(1X,A8,ES13.5))') &
    !         &   ANAME, ((TFUEL(iload)%isonam(ij), TFUEL(iload)%density(ij)),ij=kk+1,TFUEL(iload)%niso)
    !         write(UNIT,*)
    !     end do
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
            ANAME = RFUEL(iload)%name
            do ii = 1, num_pert_region
                write(UNIT,'("ADJUST_ZONE   RF_FA_",I0.5,"_",I0.2,2X,"FIRST_ORDER_PT",2X,A8,1X,A8," 1.0")') &
                &   IASS, iload, pert_regions(ii), ANAME
            end do
            write(UNIT,*)
        end do
        ! if (measure_uncertainty) then
        !     do iload = 1, NCMP2
        !         ANAME = TFUEL(iload)%name
        !         do ii = 1, num_pert_region
        !             write(UNIT,'("ADJUST_ZONE   RF_FA_",I0.5,"_",I0.2,2X,"FIRST_ORDER_PT",2X,A8,1X,A8," 1.0")') &
        !             &   IASS, iload+NCMP, pert_regions(ii), ANAME
        !         end do
        !         write(UNIT,*)
        !     end do
        ! end if
    end do
    num_pert_cases = num_pool * NCMP   ! for the time being, each assembly has at most 2 perturbations, one for fresh fuel and one for empty

    close(UNIT)

End Subroutine

