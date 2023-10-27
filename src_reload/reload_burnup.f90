Subroutine reload_burnup(NOUT)
    ! read region burnup results and compute assembly discharge burnup
    use m_reload
    use NE_Kind
    implicit none
    integer :: NOUT
    ! local
    integer :: NIN, IOS
    integer :: ireg, nreg, IZ, IX, IY, IASS
    real(8) :: total_mass
    character(len=8) :: ANAME

    100 format('[reload_burnup]...',A,1X,I6)

    NIN = NE_Kind_GetFreeLogicalUnit()
    open(NIN, file=burnup_data, status='old', form='formatted', iostat=IOS)
    if(IOS /= 0) then
        write(NOUT,100) 'Error: cannot open burnup dataset '//burnup_data
        call abort
    end if

    read(NIN,*) nreg
    region_burnup%nreg = nreg
    allocate(region_burnup%regnam(nreg), region_burnup%nburn(nreg), region_burnup%mass0(nreg), &
        region_burnup%burnup(nreg), assmas(num_xynod), assbrn(num_xynod), asscyc(num_xynod),    stat=IOS)
    if(IOS /= 0) then
        write(NOUT,100) 'Error: cannot allocate arrays for assembly burnup calculation'
        call abort
    end if
    do ireg = 1, nreg
        read(NIN,*) region_burnup%regnam(ireg), region_burnup%mass0(ireg), region_burnup%burnup(ireg), region_burnup%nburn(ireg)
    end do
    region_burnup%defined = .true.

    write(NOUT,'(/,20X,"*** Assembly burnups ***")')
    if(num_hex_sector > 0) then
        write(NOUT,'(/,"Index  IRING  IPOS   Burnup, MWD/MT    Burned Cycles",/)')
    else
        write(NOUT,'(/,"Index    IX    IY    Burnup, MWD/MT    Burned Cycles",/)')
    end if

    ! compute assembly burnups
    IASS = 0
    do IY = 1, num_ynode
        do IX = 1, num_xnode
            IASS = IASS + 1
            assbrn(IASS) = 0.0D0
            assmas(IASS) = 0.0D0
            asscyc(IASS) = 0
            if(assembly_map(IX,IY) /= reload_type) cycle    ! not fuel assembly
            do IZ = bot_active_node, top_active_node
                ireg = GEODST%MR(IX,IY,IZ)
                if(IZ > bot_active_node) THEN
                  IF (IREG == GEODST%MR(IX,IY,IZ-1)) cycle
                END IF
                !if(IZ > bot_active_node .and. ireg == GEODST%MR(IX,IY,IZ-1)) cycle  ! in case a region is over counted
                ANAME = LABELS%REGNAM(ireg)
                do ireg = 1, nreg
                    if(ANAME == region_burnup%regnam(ireg)) exit
                end do
                if(ireg > nreg) then
                    write(NOUT,100) 'warning: burnup of region '//ANAME//' is not given, will ignore it'
                else
                    assbrn(IASS) = assbrn(IASS) + region_burnup%mass0(ireg) * region_burnup%burnup(ireg)  ! assembly integrated MWD
                    assmas(IASS) = assmas(IASS) + region_burnup%mass0(ireg)
                    asscyc(IASS) = region_burnup%nburn(ireg)       ! fuel reloading is assembly based
                end if
            end do
            assbrn(IASS) = assbrn(IASS) / assmas(IASS)   ! assembly average burnup in MWD/MT
            if(num_hex_sector > 0) then
                write(NOUT,'(I5,2I6,2X,ES13.5,8X,I4)') IASS, hexagon_rp(1:2,IASS), assbrn(IASS), asscyc(IASS)
            else 
                write(NOUT,'(I5,2I6,2X,ES13.5,8X,I4)') IASS, IX, IY, assbrn(IASS), asscyc(IASS)
            end if
        end do
    end do

    ! find the average burnup of fuel assemblies
    if(min_burnup < 0.0) then
        write(NOUT,'(/,"the maximum number of assembly loading cycles =",I4)') maxval(asscyc)
        min_burnup = 0.0D0
        if(.not. manual_reload) then
            total_mass = 0.0D0
            IASS = 0
            do IY = 1, num_ynode 
                do IX = 1, num_xnode
                    IASS = IASS + 1
                    if(assembly_map(IX,IY) == reload_type) then
                        min_burnup = min_burnup + assbrn(IASS)*assmas(IASS)
                        total_mass = total_mass + assmas(IASS)
                    end if
                end do
            end do
            min_burnup = min_burnup / total_mass
            write(NOUT,'("the average burnup (MWD/MT) of fuel assemblies so far =",ES13.5)') min_burnup
        end if
    end if

End Subroutine
