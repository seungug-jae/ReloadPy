Subroutine get_boc_power(output, NOUT, keff, pwr, peak_ass, peak_pwr)
    ! read (fuel) assembly integrated powers from DIF3D output or BOC edits of REBUS output
    use m_reload
    use m_utility
    implicit none
    character(len=*) :: output        ! output file name
    integer :: NOUT, peak_ass
    real(8) :: keff, pwr(num_xynod), peak_pwr
    ! local
    integer :: NIN, IOS, nreg, ireg, idum, IX, IY, IZ, NCZ, IASS
    real(8), allocatable :: regpwr(:)
    integer, allocatable :: assreg(:)
    character(len=500) :: line
    character(len=8) :: ANAME
    real(8) :: fdum
    data NIN /100/ 
    
    100 format('[get_boc_power]...',A,1X,I6)
    open(NIN, file=output, status='old', form='formatted', iostat=IOS) 
    if(IOS /= 0) then
        write(NOUT,100) 'Error: cannot open DIF3D/REBUS output file '//output
        call abort
    end if

    ! record keff at BOC by the way
    do 
        read(NIN,'(A)',iostat=IOS) line
        if(IOS /= 0) then
            write(NOUT,100) 'Error: unexpected I/O error when searching for K-EFFECTIVE'
            call abort
        end if
        idum = index(line,'K-EFFECTIVE =')
        if(idum > 0) then
            idum = idum + 13
            read(line(idum:500),*) keff
            exit
        end if
    end do

    ! read region powers from dif3d/rebus output
    nreg = GEODST%NREG
    allocate(regpwr(nreg), assreg(GEODST%NINTK))
    ireg = 0
    loop1: do
        if(ireg == nreg) exit
        call SkipTo('0            REGION AND AREA POWER INTEGRALS FOR  K-EFF PROBLEM',.true.,NIN,IOS)
        if(IOS /= 0) then
            write(NOUT,100) 'Error: unexpected I/O error'
            call abort
        end if
        call SkipLines(NIN,4)
        do 
            read(NIN,'(A)') line
            if(len_trim(line) == 0) cycle loop1
            read(line,*) ireg, ANAME, idum, ANAME, fdum, fdum, regpwr(ireg)
            if(ireg == nreg) exit loop1
        end do
    end do loop1
    close(NIN)

    ! find peak assembly power
    peak_pwr = 0.0D0
    IASS = 0
    do IY = 1, num_ynode
        do IX = 1, num_xnode
            IASS = IASS + 1
            pwr(IASS) = 0.0D0
            if(assembly_map(IX,IY) /= reload_type) cycle
            ! locate regions in this assembly
            NCZ = 1
            assreg(1) = GEODST%MR(IX,IY,1)
            do IZ = 2, GEODST%NINTK
                if(GEODST%MR(IX,IY,IZ) == assreg(NCZ)) cycle
                NCZ = NCZ + 1
                assreg(NCZ) = GEODST%MR(IX,IY,IZ)
            end do
            ! compute assembly power
            fdum = 0.0
            do IZ = 1, NCZ
                fdum = fdum + regpwr(assreg(IZ))
            end do
            if(num_hex_sector > 0 .and. hexagon_rp(1,IASS) == 1) then  ! central hexagonal assembly
                fdum = fdum * (6.0 / num_hex_sector)
            end if
            pwr(IASS) = fdum
            ! compare peak
            if(fdum > peak_pwr) then
                peak_pwr = fdum
                peak_ass = IASS
            end if
        end do
    end do

End Subroutine