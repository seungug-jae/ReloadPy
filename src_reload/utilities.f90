Module m_utility
contains

!-----------------------------------------------------------------------
Subroutine SkipLines(NIN,N)
    ! skip the following N lines in text file NIN
    integer :: NIN, N, i, IOS
    do i = 1, N
        read(NIN,*, iostat=IOS)
        if(IOS.NE.0) then 
            write(0,'("[SkipLines]...Error: I/O error encountered when skip ",I3, &
            &   " lines following the current position of file unit",I4)') N, NIN 
            call abort 
        end if
    end do
End Subroutine

!-----------------------------------------------------------------------
Subroutine SkipTo(tag, lead, NIN, IOS)
    ! skip to the next line that contains an indicator string <tag>, 
    ! after return, I/O cursor will the at the beginning of the next line 
    implicit none
    integer :: NIN, tag_len
    character(len=*) :: tag
    logical :: lead     ! .true. - tag is the leading string in the line (ignore leading blanks); .false. - not sure
    integer :: IOS      ! < 0, end-of-file; = 0, success; > 0, I/O error
    ! local
    integer :: taglen, length
    character(len=500) :: line 
    taglen = len_trim(tag)
    
    do 
        read(NIN,'(A)',iostat=IOS) line
        if(IOS .LT. 0) then
            exit
        else if(IOS .GT. 0) then 
            write(0,'("[SkipTo]...Unexpected I/O error")')
            call abort 
        end if
        line = adjustl(line)
        length = len_trim(line)
        if(length .EQ. 0) cycle 
        if(lead) then 
            if(line(1:taglen) .EQ. tag) exit 
        else
            if(index(line(1:length),tag) .GT. 0) exit 
        end if        
    end do
    
End Subroutine

!-----------------------------------------------------------------------
Subroutine scan_arc_file_header(NIN,MCARD,NCARDS,FREEFORM,SETNAM)
    implicit none 
    integer :: NIN, MCARD, NCARDS(99), FREEFORM
    character(len=8) :: SETNAM  ! file set name (=A.NIP3, A.DIF3D, A.HMG4C, A.BURN, etc.)
    ! MCARD, NCARDS, FREEFORM are all outputs
    ! MCARD = maximum card number
    ! NCARDS = number of cards of each type
    ! FREEFORM = 1 / 0 if input is in free / fixed format
    integer :: IDUM !, ii, jj, nline

    read(NIN,'(A8,3I5/(16I9))') SETNAM, MCARD, IDUM, FREEFORM, NCARDS(1:MCARD)
    !nline = (MCARD+15)/16
    !jj = 0
    !NCARDS(:) = 0
    !do ii = 1, nline-1
    !    read(NIN,*) NCARDS(jj+1:jj+16)
    !    jj = jj + 16
    !end do
    !ii = MCARD-jj
    !read(NIN,*) NCARDS(jj+1:jj+ii)

End Subroutine

!-----------------------------------------------------------------------
! If Fortran 2008 is used, use intrinsic subroutine EXECUTE_COMMAND_LINE (command [,wait,exitstat,cmdstat,cmdmsg]) 
! instead of the following trial 
Subroutine system_command(NOUT, command)
    integer :: NOUT
    character(len=*) :: command
    ! local
    integer :: ReturnedError, I
    integer, parameter :: NonErrorSystem = -2147483648 
    integer, parameter :: NonErrorCode = 0

    ReturnedError = 0
    ReturnedError = SYSTEM(command)

    ! The SYSTEM function returns a signed integer: 
    !  1) A valid error code from the SYSTEM function is -1 (indicates that SYSTEM function itself failed)
    !  2) A non-error code from the SYSTEM function is 0, however...
    !  3) The returned integer is actually the error code from the function called by system (i.e. cp in this case) when system error=0
    ! Problems:
    ! 1) In some platforms/compilers/functions an unsigned integer rather than a signed one is returned.
    !     Clearly this is a problem if it tries to return a -1 as there is no such thing with an unsigned integer.
    !     Equivalences 0 => -2147483648 ; -1 => 4294967296 => 2147483647;
    ! 2) In some functions it is perfectly normal to recieve a > 0 value which are warnings more than errors.
    ! 3) Some implementations splice together information on a 8 bit construction of a 32 bit integer.
    !    In this case, the least important 8 bits is the returned code which for most machines

    IF (ReturnedError .EQ. NonErrorSystem)    ReturnedError = 0  ! Returned an unsigned integer value meant to be zero
    IF (ReturnedError .EQ. -NonErrorSystem-1) ReturnedError = -1 ! Tried to return -1, but unsigned integer screwed us
    IF (ReturnedError .EQ. -1) THEN
        WRITE(NOUT,'("[system_command]...Error: ",A," system command failed!",I16)') command, ReturnedError
        call abort
    ELSE IF (ReturnedError /= NonErrorCode) THEN
        !WRITE(NOUT,'("[system_command]...warning: ",A," system command returned with error code ",I16)') command, ReturnedError
        ! if going to handle problem 3), have to split 32bit integers, and turn off argument checking of compiler
        call get_8bit_errcode(ReturnedError, I)
        if(I /= 0) then
            WRITE(NOUT,'("[system_command]...Error: ",A," system command failed!")') command
            call abort
        end if
    END IF

    Contains
    Subroutine get_8bit_errcode(my8bitwords, errorcode)
        integer(kind=1) :: my8bitwords(4)
        integer :: errorcode
        errorcode = 0
        if(my8bitwords(4) == 255) errorcode = 1
    End Subroutine

End Subroutine

!-----------------------------------------------------------------------
! quick sort (unstable)
recursive subroutine quicksort(size, val, mask, first, last, order)
    implicit none
    integer :: size, first, last
    real(8) :: val(size)
    integer :: mask(size)
    character(len=*) :: order  ! 'ascend' / 'descend'
    ! local 
    real(8) :: x, t
    integer :: i, j, k

    x = val( (first+last) / 2 )
    i = first
    j = last

    if(order == 'ascend') then
        do
           do while (val(i) < x)
              i=i+1
           end do
           do while (x < val(j))
              j=j-1
           end do
           if (i >= j) exit
           t = val(i);  val(i) = val(j);  val(j) = t
           k = mask(i); mask(i) = mask(j); mask(j) = k
           i=i+1
           j=j-1
        end do
    else
        do
            do while (val(i) > x)
               i=i+1
            end do
            do while (x > val(j))
               j=j-1
            end do
            if (i >= j) exit
            t = val(i);  val(i) = val(j);  val(j) = t
            k = mask(i); mask(i) = mask(j); mask(j) = k
            i=i+1
            j=j-1
        end do
    end if
    if (first < i-1) call quicksort(size, val, mask, first, i-1, order)
    if (j+1 < last)  call quicksort(size, val, mask, j+1, last, order)

end subroutine


!-----------------------------------------------------------------------
! insertion sort (stable)
Subroutine insertSort(size, val, mask, order, tolerance)
    integer :: size
    real(8) :: val(size)
    integer :: mask(size)
    character(len=*) :: order  ! 'ascend' / 'descend'
    real(8),optional :: tolerance ! a small positive value, positive means relative tolerance, negative means absolute tolerance
    ! local
    real(8) :: tmp, tol=0.0d0
    Integer :: I, J, K

    if(present(tolerance)) tol = tolerance

    if(tol > 0.0d0) then    ! use relative difference
        if(order == 'ascend') then
            do I = 2, size
                tmp = val(I)
                K = mask(I)
                do J = I-1, 1, -1
                    If (tmp-val(J) < -tol*abs(val(J))) then
                        val(J+1) = val(J)
                        mask(J+1) = mask(J)
                    else
                        exit
                    end if
                end do
                val(J+1) = tmp
                mask(J+1) = K
            end do
        else
            do I = 2, size
                tmp = val(I)
                K = mask(I)
                do J = I-1, 1, -1
                    If (tmp-val(J) > tol*abs(val(J))) then
                        val(J+1) = val(J)
                        mask(J+1) = mask(J)
                    else
                        exit
                    end if
                end do
                val(J+1) = tmp
                mask(J+1) = K
            end do
        end if
    else    ! use absolute difference
        tol = abs(tolerance)
        if(order == 'ascend') then
            do I = 2, size
                tmp = val(I)
                K = mask(I)
                do J = I-1, 1, -1
                    If (tmp-val(J) < -tol) then
                        val(J+1) = val(J)
                        mask(J+1) = mask(J)
                    else
                        exit
                    end if
                end do
                val(J+1) = tmp
                mask(J+1) = K
            end do
        else
            do I = 2, size
                tmp = val(I)
                K = mask(I)
                do J = I-1, 1, -1
                    If (tmp-val(J) > tol) then
                        val(J+1) = val(J)
                        mask(J+1) = mask(J)
                    else
                        exit
                    end if
                end do
                val(J+1) = tmp
                mask(J+1) = K
            end do
        end if
    end if

End Subroutine


End Module