Subroutine read_boc_keff(output, NOUT, keff)
    ! read keff from DIF3D output or BOC edits of REBUS output
    implicit none
    character(len=*) :: output        ! output file name
    integer :: NOUT
    real(8) :: keff
    ! local
    integer :: NIN, IOS,idum
    character(len=500) :: line
    data NIN /100/ 
    
    100 format('[read_boc_keff]...',A,1X,I6)
    open(NIN, file=output, status='old', form='formatted', iostat=IOS) 
    if(IOS /= 0) then
        write(NOUT,100) 'Error: cannot open DIF3D/REBUS output file '//output
        call abort
    end if

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
    close(NIN)

End Subroutine