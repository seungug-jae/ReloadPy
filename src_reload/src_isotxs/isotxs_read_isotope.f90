    subroutine isotxs_read_isotope(fileunit,ngroup,user_isotope,isize_scratch,isize_scratch2, &
        scratch,scratch2,iopt)

    !  Feb. 2011 : Initial production version, C. H. Lee
    !              Modified from the UNIC routines

    !  Jan. 9 2020: modified by pdeng
    use m_kind
    use m_isotxs

    implicit none

    integer, intent(IN) :: fileunit,ngroup,iopt
    !  iopt = 0,  read record 4D and return (file is backspaced)
    !       = 1,  read record 4D and return (file is not backspaced)
    !       = 2,  start with 5D record and read all the rest data of this isotope
    !       the file I/O position is controled ouside
    type (type_isotxs_isotope), intent(INOUT) :: user_isotope
    !  scratch data
    integer, intent(IN) :: isize_scratch,isize_scratch2
    real(kind=sp), intent(OUT) :: scratch(ngroup,isize_scratch),scratch2(isize_scratch2)
    !  local
    integer :: nscmax,nsblok,currentfilerecord
    integer :: datalength
    integer :: g,i,j,k,l,m,n,nn,nl,iterm,ms,jl,ju,mcolumn,mrow,ig,ib,ie,it
    integer :: nop1,nop2,idsct_n,nwds
    logical :: Lprint_isotxs=.false.
    integer, parameter :: debug_unit = 91
    real :: sigs

    nscmax = user_isotope%nscmax
    nsblok = 1     ! no sub-blocks for scattering matrix any more
    currentfilerecord = 0

    if(iopt < 2) then
        !     card type 4: isotope information
        read(fileunit)  user_isotope%habsid,user_isotope%hident,user_isotope%hmat,    &
            user_isotope%amass,user_isotope%efiss,user_isotope%ecapt,     &
            user_isotope%temp,user_isotope%sigpot,user_isotope%adens,     &
            user_isotope%kbr,user_isotope%ichi,user_isotope%ifis,         &
            user_isotope%ialf,user_isotope%inp,user_isotope%in2n,         &
            user_isotope%ind,user_isotope%int,user_isotope%ltot,          &
            user_isotope%ltrn,user_isotope%istrpd,                        &
            (user_isotope%idsct(n),n=1,nscmax),                           &
            (user_isotope%lord(n),n=1,nscmax),                            &
            ((user_isotope%jband(g,n),g=1,ngroup),n=1,nscmax),            &
            ((user_isotope%ijj(g,n),g=1,ngroup),n=1,nscmax)

        nwds=1+user_isotope%ltrn+user_isotope%ltot+user_isotope%ialf+user_isotope%inp &
            +user_isotope%in2n+user_isotope%istrpd+2*user_isotope%ifis &
            +user_isotope%ichi*(2/(user_isotope%ichi+1))
        !     if(nwds > (10+2*user_isotope%maxord)) then
        !       write(*,*) 'nwds > (10+2*maxord)=',nwds,(10+2*user_isotope%maxord)
        !     endif
        if(iopt == 0) backspace fileunit
        return
    end if

    currentfilerecord = currentfilerecord + 1 ! increment the record position

    !  card type 5: principal cross sections
    !  nop is the number of moments in the cross section data, capture is always present
    nop1 = isotxs_isotope_calculatenop(user_isotope,50) ! this will add the principle moments from 1 to 11
    read(fileunit)  ((scratch(g,i),g=1,ngroup),i=1,nop1)
    currentfilerecord = currentfilerecord + 1

    !  card type 5: principal cross sections
    !  transport xs data
    nop1 = 0
    nop2 = isotxs_isotope_calculatenop(user_isotope,1)
    if(nop2 > nop1) then
        do i = 1,(nop2-nop1)
            user_isotope%transport(1:ngroup,i) = scratch(1:ngroup,nop1+i)
        enddo
    endif

    !  total xs data
    nop1 = nop2
    nop2 = isotxs_isotope_calculatenop(user_isotope,2)
    if(nop2 > nop1) then
        do i = 1,(nop2-nop1)
            user_isotope%total(1:ngroup,i) = scratch(1:ngroup,nop1+i)
        enddo
    endif

    !  n,gamma xs data
    nop1 = nop2
    nop2 = isotxs_isotope_calculatenop(user_isotope,3)
    if(nop2 > nop1) user_isotope%ngamma(1:ngroup) = scratch(1:ngroup,nop2)

    !  n,fission xs data
    nop1 = nop2
    nop2 = isotxs_isotope_calculatenop(user_isotope,4)
    if(nop2 > nop1) user_isotope%fission(1:ngroup) = scratch(1:ngroup,nop2)

    !  nu xs data
    nop1 = nop2
    nop2 = isotxs_isotope_calculatenop(user_isotope,5)
    if(nop2 > nop1) user_isotope%nu(1:ngroup) = scratch(1:ngroup,nop2)

    !  chi xs data
    nop1 = nop2
    nop2 = isotxs_isotope_calculatenop(user_isotope,6)
    if(nop2 > nop1) user_isotope%chi(1:ngroup,1) = scratch(1:ngroup,nop2)

    !  n,alpha xs data
    nop1 = nop2
    nop2 = isotxs_isotope_calculatenop(user_isotope,7)
    if(nop2 > nop1) user_isotope%alpha(1:ngroup) = scratch(1:ngroup,nop2)

    !  n,proton xs data
    nop1 = nop2
    nop2 = isotxs_isotope_calculatenop(user_isotope,8)
    if(nop2 > nop1) user_isotope%proton(1:ngroup) = scratch(1:ngroup,nop2)

    !  n,2n xs data
    nop1 = nop2
    nop2 = isotxs_isotope_calculatenop(user_isotope,9)
    if(nop2 > nop1) user_isotope%n2n(1:ngroup) = scratch(1:ngroup,nop2)

    !  n,deuteron xs data
    nop1 = nop2
    nop2 = isotxs_isotope_calculatenop(user_isotope,10)
    if(nop2 > nop1) user_isotope%deuteron(1:ngroup) = scratch(1:ngroup,nop2)

    !  n,tritium xs data
    nop1 = nop2
    nop2 = isotxs_isotope_calculatenop(user_isotope,11)
    if(nop2 > nop1) user_isotope%tritium(1:ngroup) = scratch(1:ngroup,nop2)

    !  coordinate directed xs data
    nop1 = nop2
    nop2 = isotxs_isotope_calculatenop(user_isotope,12)
    if(nop2 > nop1) then
        do i = 1,(nop2-nop1)
            user_isotope%strpd(1:ngroup,i) = scratch(1:ngroup,nop1+1)
        enddo
    endif

    if(Lprint_isotxs) then
        write(*,*) 'printing xs ... ', user_isotope%habsid, user_isotope%hisonm
        write(debug_unit,*) user_isotope%habsid
        write(debug_unit,'(5x,3a12,$)') 'sigtr','sigt','sigc'
        if(user_isotope%ifis > 0) write(debug_unit,'(a12,$)') 'sigf','nu','chi'
        if(user_isotope%in2n > 0) write(debug_unit,'(a12,$)') 'sign2n'
        if(user_isotope%ialf > 0) write(debug_unit,'(a12,$)') 'signalpha'
        if(user_isotope%inp  > 0) write(debug_unit,'(a12,$)') 'signp'
        if(user_isotope%ind  > 0) write(debug_unit,'(a12,$)') 'signd'
        if(user_isotope%int  > 0) write(debug_unit,'(a12,$)') 'signt'
        write(debug_unit,*)
        do ig=1,ngroup
            write(debug_unit,'(i5,1p3e12.4,$)') ig,user_isotope%transport(ig,1),user_isotope%total(ig,1),user_isotope%ngamma(ig)
            if(user_isotope%ifis > 0) write(debug_unit,'(1p3e12.4,$)') user_isotope%fission(ig),user_isotope%nu(ig),user_isotope%chi(ig,1)
            if(user_isotope%in2n > 0) write(debug_unit,'(1pe12.4,$)')  user_isotope%n2n(ig)
            if(user_isotope%ialf > 0) write(debug_unit,'(1pe12.4,$)')  user_isotope%alpha(ig)
            if(user_isotope%inp  > 0) write(debug_unit,'(1pe12.4,$)')  user_isotope%proton(ig)
            if(user_isotope%ind  > 0) write(debug_unit,'(1pe12.4,$)')  user_isotope%deuteron(ig)
            if(user_isotope%int  > 0) write(debug_unit,'(1pe12.4,$)')  user_isotope%tritium(ig)
            write(debug_unit,*)
        enddo
        if(user_isotope%ltot > 1 .or. user_isotope%ltrn > 1) then
            write(debug_unit,'(5x,2a20)') '(sigt L=2,N)','(sigtr L=2,N-1)'
            do ig=1,ngroup        ! user_isotope%ltrn = user_isotope%ltot - 1
                write(debug_unit,'(i5,1p100e12.4)') ig,(user_isotope%total(ig,L),L=2,user_isotope%ltot),(user_isotope%transport(ig,L),L=2,user_isotope%ltrn)
            enddo
        endif
    endif

    !  card type 6: isotope chi data
    if(user_isotope%ichi > 1) then
        nop1 = user_isotope%ichi
        read(fileunit)  ((user_isotope%chi(g,i),i=1,nop1),g=1,ngroup),(user_isotope%isopec(g),g=1,ngroup)
        currentfilerecord = currentfilerecord + 1
    endif

    !  card type 7: scattering sub-bloc        
    do n=1,nscmax                                               ! nscmax is the total number of blocked data, each block can have subblocks
        if(user_isotope%lord(n) > 0) then                        ! lord(n) = number of subblocks in block, (0 = none, 1 = p0, 2 = p1, ...)
            idsct_n=user_isotope%idsct(n)/100                     ! idsct lists the block type; idsct_n = 0,1,2,3,4,...
            ! 000 is total, 100 is elastic, 200 is inelastic, 300 is n2n, all else is ignored
            nn=user_isotope%idsct(n)-idsct_n*100                  ! nn is the starting anisotropy order of the subblock; 0,1,2,3,4... p0, p1, p2, ...
            do ms=1,nsblok                                        ! ms = 1 always, nsblock > 1 is used for partitioning of subblocks
                jl=(ms-1)*((ngroup-1)/nsblok+1)+1                  ! jl = 1 always
                ju=ms*((ngroup-1)/nsblok+1)                        ! ju = ngroup always
                datalength = 0
                do j=jl,ju                                         ! j=1,ngroup
                    datalength=datalength+user_isotope%jband(j,n)    ! datalength = length of single subblock data (all subblocks have the same form)
                enddo
                m = datalength
                !           write(iunit,*)'scat data length = ',m
                datalength=datalength*user_isotope%lord(n)         ! datalength = number of subblocks * length of subblock = all data in block
                if(datalength /= 0) then
                    read(fileunit)  scratch2(1:datalength)
                    !              now seperate the data and dissemenate into the data arrays
                    l = 0                                           ! data index of all subblocked data
                    do nl=1,user_isotope%lord(n)                    ! nl = 1,number of subblocks; extract each subblock in this block
                        iterm = 0
                        if(idsct_n == 0) iterm = user_isotope%total_start
                        if(idsct_n == 1) iterm = user_isotope%elast_start
                        if(idsct_n == 2) iterm = user_isotope%inela_start
                        if(idsct_n == 3) iterm = user_isotope%n2n_start
                        if(iterm > 0) iterm = iterm + nn + nl - 1    ! starting location for iterm

                        !                 assuming a lower triangular setup for the scattering data in scattering(g,g') where g'-> g
                        !                 the following data parses that array into rows which are stored contiguously from row 1 to row g, the right most column to the left most column
                        do j=jl,ju
                            mcolumn = j                                ! defines the current group crap is scattering into
                            if(user_isotope%jband(j,n) > 0) then       ! if = 0 then no data exists
                                do k = 1,user_isotope%jband(j,n)        ! transfer this piece of the subblocked data
                                    l = l + 1                            ! increment the index of the read in data
                                    mrow = (j+user_isotope%ijj(j,n)) - k ! bottom to top assumption: gives the group from which scattering is coming from
                                    if(iterm > 0) user_isotope%scattering(j,mrow,iterm) = scratch2(l)
                                enddo
                            endif
                        enddo  ! j=jl,ju
                    enddo  ! nl=1,user_isotope%lord(n)
                    currentfilerecord = currentfilerecord + 1
                endif  ! if(datalength /= 0)
            enddo  ! ms=1,nsblok
        endif  ! if(user_isotope%lord(n) > 0)

        if(Lprint_isotxs) then
            if(idsct_n == 0) then
                write(debug_unit,*) 'total scattering matrix',user_isotope%idsct(n),mod(user_isotope%idsct(n),100)
                do mrow=1,ngroup
                    sigs=0.
                    do j=1,ngroup
                        sigs=sigs+user_isotope%scattering(j,mrow,iterm)
                    enddo
                    call groupbound(ngroup,user_isotope%scattering(:,mrow,iterm),ib,ie)
                    it=max(ie-ib+1,0)
                    write(debug_unit,'(4i5,1p3000e12.4)') mrow,ib,ie,it,sigs,(user_isotope%scattering(j,mrow,iterm),j=ib,ie)
                enddo
            endif
            if(idsct_n == 1) then
                write(debug_unit,*) 'elastic scattering matrix',user_isotope%idsct(n),mod(user_isotope%idsct(n),100)
                do mrow=1,ngroup
                    sigs=0.
                    do j=1,ngroup
                        sigs=sigs+user_isotope%scattering(j,mrow,iterm)
                    enddo
                    call groupbound(ngroup,user_isotope%scattering(:,mrow,iterm),ib,ie)
                    it=max(ie-ib+1,0)
                    write(debug_unit,'(4i5,1p3000e12.4)') mrow,ib,ie,it,sigs,(user_isotope%scattering(j,mrow,iterm),j=ib,ie)
                enddo
            endif
            if(idsct_n == 2) then
                write(debug_unit,*) 'inelastic scattering matrix',user_isotope%idsct(n),mod(user_isotope%idsct(n),100)
                do mrow=1,ngroup
                    sigs=0.
                    do j=1,ngroup
                        sigs=sigs+user_isotope%scattering(j,mrow,iterm)
                    enddo
                    call groupbound(ngroup,user_isotope%scattering(:,mrow,iterm),ib,ie)
                    it=max(ie-ib+1,0)
                    write(debug_unit,'(4i5,1p3000e12.4)') mrow,ib,ie,it,sigs,(user_isotope%scattering(j,mrow,iterm),j=ib,ie)
                enddo
            endif
            if(idsct_n == 3) then
                write(debug_unit,*) 'n2n scattering matrix',user_isotope%idsct(n),mod(user_isotope%idsct(n),100)
                do mrow=1,ngroup
                    sigs=0.
                    do j=1,ngroup
                        sigs=sigs+user_isotope%scattering(j,mrow,iterm)
                    enddo
                    call groupbound(ngroup,user_isotope%scattering(:,mrow,iterm),ib,ie)
                    it=max(ie-ib+1,0)
                    write(debug_unit,'(4i5,1p3000e12.4)') mrow,ib,ie,it,sigs,(user_isotope%scattering(j,mrow,iterm),j=ib,ie)
                enddo
            endif
        endif

    enddo  ! n=1,nscmax

    end subroutine



    subroutine groupbound(ng,sig,ib,ie)

    implicit none

    integer, intent(in) :: ng
    integer, intent(out) :: ib,ie
    real, intent(in) :: sig(ng)

    !  Local
    integer j

    do j=1,ng
        if(sig(j) > 0.) exit
    enddo
    ib=j
    do j=ng,1,-1
        if(sig(j) > 0.) exit
    enddo
    ie=j

    end

