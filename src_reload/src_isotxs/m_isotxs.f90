    module m_isotxs

    !  Feb. 2011 : Initial production version, C. H. Lee
    !              Modified from the UNIC routines
    !  May 2021: Modified by P. Deng
    !
    ! === module functions:
    !  function   isotxs_isotope_calculatenop(user_isotope,selectedmoment)
    !  subroutine isotxs_allocate_header(niso,ng,ichist,user_isotxs)
    !  subroutine isotxs_deallocate_header(user_isotxs)
    !  subroutine isotxs_allocate_isotope(user_isotope,ngroup,nscmax,user_lord,user_idsct)
    !  subroutine isotxs_deallocate_isotope(user_isotope)
    !  subroutine isotxs_isotope_computeband(user_isotope)
    ! ===> new routines
    !  subroutine isotxs_calc_group_vel(ngroup, energy, vel)
    !  subroutine isotxs_isotope_computeMAXORD(User_Isotope)
    !  subroutine isotxs_isotope_computeMAXUP(User_Isotope, MAXUP, MAXDN)
    !  function   isotxs_isotope_countrecords(isotope) result(num_records)
    !  subroutine isotxs_xs_scale_uniform(isotope, weight, newiso)
    !  subroutine isotxs_xs_scale_groupwise(isotope, ngroup, lord, weight, newiso)
    !  subroutine isotxs_get_total_scatmx(isotope, ngroup, maxord, scat)
    !  subroutine isotxs_convert_scatmx_to_total(isotxs, isotope)   -> incomplete
    !  function   isotxs_isotope_n2nfactor(isotope) result(real::n2nfactor)
    !  subroutine isotxs_transport_scatgg(isotope)
    !  subroutine isotxs_revert_transport_scatgg(isotope)


    use m_kind

    implicit none

    !  for card type 4,5,6,7
    type type_isotxs_isotope
        !    data information needed for allocation or copy of file header
        logical  :: defined=.false.                       ! indicates that isotope data is fully defined (type 4,5,6 and 7)
        logical  :: t04_defined=.false.                   ! indicates that card type 4 data was allocated
        integer  :: ngroup = 0                            ! number of groups
        integer  :: nscmax = 0                            ! maximum number of blocks of scattering data per isotope
        integer  :: maxord = 0                            ! maximum legendre order for this isotope
        integer  :: maxterms=0                            ! maximum number of scattering terms = (sum(lord(n)),n=1,nscmax). Given lord(n) = 1 for each block, maxterms = nscmax
        character(len=8) :: hisonm = ''                   ! copy of file wide data
        logical  :: scatgg_corrected = .false.            ! pdeng: whether extended transport approximation is made to self-scattering (diagonal terms)
        !    card type 4: isotope control and group independent data
        character(len=8) :: habsid = ''                   ! hollerith absolute isotope label - same for all versions of the same isotope in file (a6)
        character(len=8) :: hident = ''                   ! identifier of library from which basic data came(e.g. endf/b) (a6)
        character(len=8) :: hmat   = ''                   ! isotope identification (e.g. endf/b mat no.) (a6)
        real(kind=sp) :: amass  = 0.0                     ! gram atomic weight
        real(kind=sp) :: efiss  = 0.0                     ! total thermal energy yield/fission (w.sec/fiss)
        real(kind=sp) :: ecapt  = 0.0                     ! total thermal energy yield/capture (w.sec/capt)
        real(kind=sp) :: temp   = 0.0                     ! isotope temperature (degrees kelvin)
        real(kind=sp) :: sigpot = 0.0                     ! average effective potential scattering in resonance range (barns/atom)
        real(kind=sp) :: adens  = 0.0                     ! density of isotope in mixture in which isotope cross sections were generated (a/barn-cm)
        integer :: kbr    = 0                             ! isotope classification: 0=undefined        1=fissile    2=fertile     3=other actinide
        !                         4=fission product  5=structure  6=coolant     7=control
        integer :: ichi   = 0                             ! isotope fission spectrum flag
        ! ichi=0    use file-wide chi
        ! ichi=1    isotope chi vector
        ! ichi>1    isotope chi matrix
        integer :: ifis   = 0                             ! (n,f)          flag (0/1)=data (isn't,is) present
        integer :: ialf   = 0                             ! (n,alpha)      flag (0/1)=data (isn't,is) present
        integer :: inp    = 0                             ! (n,proton)     flag (0/1)=data (isn't,is) present
        integer :: in2n   = 0                             ! (n,2n)         flag (0/1)=data (isn't,is) present
        integer :: ind    = 0                             ! (n,deuteron)   flag (0/1)=data (isn't,is) present
        integer :: int    = 0                             ! (n,tritium)    flag (0/1)=data (isn't,is) present
        integer :: ltot   = 0                             ! (n,total)      moments (0/n)=data (isn't present,has n moments)
        integer :: ltrn   = 0                             ! (n,transport)  moments (0/n)=data (isn't present,has n moments)
        integer :: istrpd = 0                             ! (n,directionl) moments (0/n)=data (isn't present,has n moments)
        integer :: iine   = 0                             ! inelastic scat flag (0/1)=data (isn't,is) present
        integer, pointer :: idsct(:)                      ! (nscmax)
        ! scattering matrix type identification for scattering block n
        ! idsct(n)=000 + nn, total scattering, (sum of elastic,inelastic, and n,2n scattering matrix terms)
        !         =100 + nn, elastic scattering
        !         =200 + nn, inelastic scattering
        !         =300 + nn, (n,2n) scattering,----see note below----
        ! where nn is the legendre expansion index of the first matrix in block n
        integer, pointer :: lord(:)                       ! (nscmax)
        ! number of scattering orders in block n.
        ! if lord(n)=0, this block is not present for this isotope
        ! the matrices in this block have legendre expansion indices of nn,nn+1,nn+2,...,nn+lord(n)-1
        integer, pointer :: jband(:,:)                    ! (ngroup,nscmax)
        ! number of groups that scatter into group j, including self-scatter, in scattering block n.
        !                 if jband(j,n)=0, no scatter data is present in block n
        integer, pointer :: ijj(:,:)                      ! (ngroup,nscmax)
        ! position of in-group scattering cross section in scattering data for group j,
        ! scattering block n,counted from the first word of group j data.
        ! if jband(j,n)/=0 then ijj(j,n) must satisfy the relation 1 <= ijj(j,n) <= jband(j,n)
        !    card type 5: principal cross sections
        real(kind=sp),pointer :: transport(:,:)           ! (ngroup,ltrn)  p(l) weighted transport cross section.
        !                      the first vector is the current (p1) weighted transport cross section.
        !                      the legendre expansion coefficient factor (2l+1) is not included.
        real(kind=sp),pointer :: total(:,:)               ! (ngroup,ltot)  p(l) weighted total cross section.
        !                      the first vector is the flux (p0) weighted total cross section.
        !                      the legendre expansion coefficient factor (2l+1) is not included.
        real(kind=sp),pointer :: ngamma(:)                ! (ngroup)       (n,gamma) cross section
        real(kind=sp),pointer :: fission(:)               ! (ngroup)       (n,fission) cross section
        real(kind=sp),pointer :: nu(:)                    ! (ngroup)       nu- number of neutrons produced per fission per group g
        real(kind=sp),pointer :: chi(:,:)                 ! (ngroup,ichi)  chi- fraction of neutronxs emitted into group g from fission in
        real(kind=sp),pointer :: alpha(:)                 ! (ngroup)       (n,alpha) cross section
        real(kind=sp),pointer :: proton(:)                ! (ngroup)       (n,proton) cross section
        real(kind=sp),pointer :: n2n(:)                   ! (ngroup)       (n,2n) cross section. where = 0.5*(sum of scat(j to g) summed over all g)
        !                      flux in group j gives the rate at which n,2n reactions occur in group j.
        real(kind=sp),pointer :: deuteron(:)              ! (ngroup)       (n,deuteron)
        real(kind=sp),pointer :: tritium(:)               ! (ngroup)       (n,tritium)
        real(kind=sp),pointer :: strpd(:,:)               ! (ngroup,istrpd) coordinate direction moments of the transport cross section
        !    the scattering data is stored in full block form rather than decomposed form. the translation is still maintained however for reading and writing.
        !    card type 6: isotope chi data
        !    data is stored automatically in chi
        integer, pointer :: isopec(:)                     ! (ngroup)   isopec(i)=k implies that spectrum k is used to calculate emission spectrum from fission in group i

        !    card type 7: scattering sub-block
        !    idsct lists the block type; idsct(n) = 0,1,2,3,4,...nscmax   000 is total, 100 is elastic, 200 is inelastic, 300 is n2n
        !    lord gives the number of anisotropic order blocks in current block of scattering data
        integer :: total_start=0, total_pnterms=0         ! starting location/number of pn terms in scattering of total cross section data
        integer :: elast_start=0, elast_pnterms=0         ! starting location/number of pn terms in scattering of elastic cross section data
        integer :: inela_start=0, inela_pnterms=0         ! starting location/number of pn terms in scattering of inelastic cross section data
        integer :: n2n_start=0,   n2n_pnterms=0           ! starting location/number of pn terms in scattering of n2n cross section data
        real(kind=sp), pointer :: scattering(:,:,:)       ! ngroup,ngroup,maxterms  (g,g',l) where g'-> g is the assumed rule
    end type

    !  this data type contains all of the data in full block form
    type type_isotxs
        logical    :: defined=.false.                     ! indicates that data structure is defined
        logical    :: n2nfactorchecked = .false.          ! indicates that the n2n factor is known for this data type
        logical    :: load_all_isotopes = .true.          ! indicates whether all isotope data are loaded into memory
        integer    :: file_unit				              ! used only when load_all_isotopes = .false.
        integer    :: nskip_header = 0                    ! used only when load_all_isotopes = .false.
        real(kind=sp) :: n2nfactor = 1.0                  ! the n2n factor (either a 1.0 or 2.0)
        !    card type 0: file identification
        character(len=8) :: hname   = 'ISOTXS'
        character(len=8) :: huse(2) = ''
        integer  :: ivers     = 0
        !    card type 1: file control
        integer  :: ngroup  = 0                           ! number of energy groups in file
        integer  :: niso    = 0                           ! number of isotopes in file
        integer  :: maxup   = 0                           ! maximum number of upscatter groups
        integer  :: maxdn   = 0                           ! maximum number of downscatter groups
        integer  :: maxord  = 0                           ! maximum scattering order (maximum value of legendre expansion index used in file)
        integer  :: ichist  = 0                           ! file-wide fission spectrum flag
        ! ichist  = 0  ! no file-wide spectrum
        ! ichist  = 1  ! file-wide chi vector
        ! ichist  > 1  ! file-wide chi matrix
        integer  :: nscmax  = 0                           ! maximum number of blocks of scattering data
        integer  :: nsblok  = 0                           ! subblocking control for scatter matrices.
        ! in the olden days the scattering data are subblocked into nsblok records(subblocks) per scattering block.
        !    card type 2: file data
        character(len=8) :: hsetid(12) = ''               ! hollerith identification of file (a6)
        character(len=8), pointer :: hisonm(:)            ! (niso)          hollerith isotope label for isotope i (a6)
        real(kind=sp), pointer :: filechi(:,:)            ! (ngroup,ichist) file-wide fission spectrum(present if ichist==1)
        real(kind=sp), pointer :: vel(:)                  ! (ngroup)        mean neutron velocity in group j (cm/sec)
        real(kind=sp), pointer :: energy(:)               ! (ngroup+1)      energy group bounds for groups 1 down to g last number is bottom
        integer, pointer :: loca(:)                       ! (niso)          number of records to be skipped to read data for isotope i
        !                 the first number is 0 indicating that the first isotope is at record position 0
        !                 note that you could fuck with this to place isotopes anywhere you want on the file
        !    card type 3: file-wide chi data
        !    chi data is stored in filechi
        integer, pointer :: isspec(:)                     ! (ngroup)        isspec(i)=k implies that spectrum k is used to calculate emission spectrum from fission in group i

        !    card type 4,5,6,7
        type (type_isotxs_isotope), pointer :: isotopedata(:) ! (niso) isotope based cross section data
    end type

    interface isotxs_xs_scale    ! generic
    module procedure isotxs_xs_scale_groupwise, isotxs_xs_scale_uniform
    end interface isotxs_xs_scale

    contains
    
    !---------------------------------------------------------------------------------------------------
    !  integer function isotxs_isotope_calculatenop(user_isotope,selectedmoment)
    !  function that calculates the number of principle moments in the data structure
    !---------------------------------------------------------------------------------------------------
    integer function isotxs_isotope_calculatenop(user_isotope,selectedmoment)

    type (type_isotxs_isotope), intent(IN) :: user_isotope
    integer, intent(IN) :: selectedmoment
    !  local
    integer nop

    nop = 0
    if(selectedmoment >=  1) nop = nop + user_isotope%ltrn
    if(selectedmoment >=  2) nop = nop + user_isotope%ltot
    if(selectedmoment >=  3) nop = nop + 1                 ! (n,gamma)
    if(selectedmoment >=  4) nop = nop + user_isotope%ifis ! fission
    if(selectedmoment >=  5) nop = nop + user_isotope%ifis ! nu
    if(selectedmoment >=  6) nop = nop + user_isotope%ichi*(2/(user_isotope%ichi+1))
    if(selectedmoment >=  7) nop = nop + user_isotope%ialf
    if(selectedmoment >=  8) nop = nop + user_isotope%inp
    if(selectedmoment >=  9) nop = nop + user_isotope%in2n
    if(selectedmoment >= 10) nop = nop + user_isotope%ind
    if(selectedmoment >= 11) nop = nop + user_isotope%int
    if(selectedmoment >= 12) nop = nop + user_isotope%istrpd

    isotxs_isotope_calculatenop = nop

    end function

    !---------------------------------------------------------------------------------------------------
    subroutine isotxs_allocate_header(niso,ng,ichist,user_isotxs)

    integer, intent(IN) :: niso,ng,ichist
    type(type_isotxs), intent(INOUT) :: user_isotxs

    !  local
    integer :: ig

    user_isotxs%niso  =niso
    user_isotxs%ngroup=ng
    user_isotxs%ichist=ichist

    allocate(user_isotxs%hisonm(niso),user_isotxs%vel(ng),user_isotxs%energy(ng+1), &
        user_isotxs%loca(niso))
    do ig=1,ng
        user_isotxs%vel(ig)=0.
    enddo

    if(ichist > 0) then
        allocate(user_isotxs%filechi(ng,ichist))
    endif

    end subroutine

    !---------------------------------------------------------------------------------------------------
    subroutine isotxs_deallocate_header(user_isotxs)

    type(type_isotxs) :: user_isotxs

    deallocate(user_isotxs%hisonm, user_isotxs%vel, user_isotxs%energy, user_isotxs%loca, user_isotxs%isotopedata)
    if(associated(user_isotxs%filechi)) then
        deallocate(user_isotxs%filechi)
    endif

    end subroutine

    !---------------------------------------------------------------------------------------------------
    subroutine isotxs_allocate_isotope(user_isotope,ngroup,nscmax,user_lord,user_idsct)
    integer, intent(IN) :: ngroup,nscmax
    integer, intent(IN) :: user_lord(nscmax),user_idsct(nscmax)       ! this is lord for each scattering block
    type (type_isotxs_isotope), intent(INOUT) :: user_isotope         ! a user data variable to be defined by reading in a isotxs file
    !  local
    integer ios,n,ilord,idsct_n,nn

    !  data for card type 4
    allocate(user_isotope%idsct(nscmax),user_isotope%lord(nscmax),                     &
        user_isotope%jband(ngroup,nscmax),user_isotope%ijj(ngroup,nscmax),stat=ios)
    user_isotope%ngroup   = ngroup
    user_isotope%nscmax   = nscmax
    if(nscmax > 0) then
        user_isotope%idsct = 0
        user_isotope%lord  = 0
        user_isotope%jband = 0
        user_isotope%ijj   = 0
    endif
    !  data for card type 5 ~ 7
    user_isotope%lord(1:nscmax) = user_lord(1:nscmax)
    user_isotope%idsct(1:nscmax) = user_idsct(1:nscmax)
    call isotxs_allocate_isotope_xs(user_isotope, ngroup, nscmax, user_lord, user_idsct)
    user_isotope%defined  = .true.

    end subroutine


    !---------------------------------------------------------------------------------------------------
    subroutine isotxs_allocate_isotope_xs(user_isotope,ngroup,nscmax,user_lord,user_idsct)
    integer, intent(IN) :: ngroup,nscmax
    integer, intent(IN) :: user_lord(nscmax),user_idsct(nscmax)       ! this is lord for each scattering block
    type (type_isotxs_isotope), intent(INOUT) :: user_isotope         ! a user data variable to be defined by reading in a isotxs file
    !  local
    integer ios,n,ilord,idsct_n,nn

    !  ngamma data
    allocate(user_isotope%ngamma(user_isotope%ngroup),stat=ios)
    user_isotope%ngamma(:) = 0.

    !  chi data
    if(user_isotope%ichi == 1) then
        allocate(user_isotope%chi(user_isotope%ngroup,user_isotope%ichi),stat=ios)
        user_isotope%chi(:,1) = 0.
    else if(user_isotope%ichi > 1) then
        allocate(user_isotope%chi(user_isotope%ngroup,user_isotope%ichi),user_isotope%isopec(user_isotope%ngroup),stat=ios)
        user_isotope%chi(:,:) = 0.
        user_isotope%isopec(:) = 0
    endif

    !  fission & nu data
    if(user_isotope%ifis > 0) then
        allocate(user_isotope%fission(user_isotope%ngroup),user_isotope%nu(user_isotope%ngroup),stat=ios)
        user_isotope%fission(:) = 0.
        user_isotope%nu(:) = 0.
    endif

    !  n,alpha data
    if(user_isotope%ialf > 0) then
        allocate(user_isotope%alpha(user_isotope%ngroup),stat=ios)
        user_isotope%alpha(:) = 0.
    endif

    !  n,proton data
    if(user_isotope%inp > 0) then
        allocate(user_isotope%proton(user_isotope%ngroup),stat=ios)
        user_isotope%proton(:) = 0.
    endif

    !  n,2n data
    if(user_isotope%in2n > 0) then
        allocate(user_isotope%n2n(user_isotope%ngroup),stat=ios)
        user_isotope%n2n(:) = 0.
    endif

    !  n,deuteron data
    if(user_isotope%ind > 0) then
        allocate(user_isotope%deuteron(user_isotope%ngroup),stat=ios)
        user_isotope%deuteron(:) = 0.
    endif

    !  n,tritium data
    if(user_isotope%int > 0) then
        allocate(user_isotope%tritium(user_isotope%ngroup),stat=ios)
        user_isotope%tritium(:) = 0.
    endif

    !  transport data
    if(user_isotope%ltrn > 0) then
        allocate(user_isotope%transport(user_isotope%ngroup,user_isotope%ltrn),stat=ios)
        user_isotope%transport(:,:) = 0.0
    endif

    !  total data
    if(user_isotope%ltot > 0) then
        allocate(user_isotope%total(user_isotope%ngroup,user_isotope%ltot),stat=ios)
        user_isotope%total(:,:) = 0.0
    endif

    !  strpd data
    if(user_isotope%istrpd > 0) then
        allocate(user_isotope%strpd(user_isotope%ngroup,user_isotope%istrpd),stat=ios)
        user_isotope%strpd(:,:) = 0.0
    endif

    !  scattering data
    !  before allocating the scattering we must determine the maximum number of terms needed to store it
    user_isotope%total_start = 0
    user_isotope%elast_start = 0
    user_isotope%inela_start = 0
    user_isotope%n2n_start   = 0
    user_isotope%total_pnterms = 0
    user_isotope%elast_pnterms = 0
    user_isotope%inela_pnterms = 0
    user_isotope%n2n_pnterms   = 0

    do n = 1,user_isotope%nscmax
        ilord = user_isotope%lord(n)
        if(ilord > 0) then         ! lord(n) is the number of subblocks of data
            idsct_n=user_isotope%idsct(n)/100         ! idsct lists the block type; idsct_n = 0,1,2,3,4,...
            nn=user_isotope%idsct(n)-idsct_n*100      ! nn is the starting anisotropy order of the subblock; 0,1,2,3,4... p0, p1, p2, ...
            if((idsct_n == 0) .and. (ilord+nn > user_isotope%total_pnterms)) user_isotope%total_pnterms = ilord + nn
            if((idsct_n == 1) .and. (ilord+nn > user_isotope%elast_pnterms)) user_isotope%elast_pnterms = ilord + nn
            if((idsct_n == 2) .and. (ilord+nn > user_isotope%inela_pnterms)) user_isotope%inela_pnterms = ilord + nn
            if((idsct_n == 3) .and. (ilord+nn > user_isotope%n2n_pnterms)  ) user_isotope%n2n_pnterms   = ilord + nn
        endif
    enddo

    user_isotope%maxterms = 0
    if(user_isotope%total_pnterms > 0) user_isotope%total_start = 1
    user_isotope%maxterms = user_isotope%maxterms + user_isotope%total_pnterms
    if(user_isotope%elast_pnterms > 0) user_isotope%elast_start = user_isotope%maxterms + 1
    user_isotope%maxterms = user_isotope%maxterms + user_isotope%elast_pnterms
    if(user_isotope%inela_pnterms > 0) user_isotope%inela_start = user_isotope%maxterms + 1
    user_isotope%maxterms = user_isotope%maxterms + user_isotope%inela_pnterms
    if(user_isotope%n2n_pnterms   > 0) user_isotope%n2n_start   = user_isotope%maxterms + 1
    user_isotope%maxterms = user_isotope%maxterms + user_isotope%n2n_pnterms

    if(user_isotope%maxterms > 0) then
        allocate(user_isotope%scattering(user_isotope%ngroup,user_isotope%ngroup,user_isotope%maxterms),stat=ios)
        if(ios /= 0) then
            write(0,'("[isotxs_allocate_isotope_xs]...error: failed to allocate storage for scattering matrices of isotope ",A)') user_isotope%hisonm
            call abort('')
        end if
        user_isotope%scattering(:,:,:) = 0.
    endif

    end subroutine

    !---------------------------------------------------------------------------------------------------
    subroutine isotxs_deallocate_isotope(user_isotope)

    !  a user data variable to be defined by reading in a isotxs file
    type (type_isotxs_isotope), intent(INOUT) :: user_isotope

    deallocate(user_isotope%idsct,user_isotope%lord,user_isotope%jband,user_isotope%ijj)

    !  ngamma data
    deallocate(user_isotope%ngamma)

    !  chi data
    if(user_isotope%ichi == 1) then
        deallocate(user_isotope%chi)
    else if(user_isotope%ichi > 1) then
        deallocate(user_isotope%chi,user_isotope%isopec)
    endif

    !  fission & nu data
    if(user_isotope%ifis > 0) deallocate(user_isotope%fission,user_isotope%nu)

    !  n,alpha data
    if(user_isotope%ialf > 0) deallocate(user_isotope%alpha)

    !  n,proton data
    if(user_isotope%inp > 0) deallocate(user_isotope%proton)

    !  n,2n data
    if(user_isotope%in2n > 0) deallocate(user_isotope%n2n)

    !  n,deuteron data
    if(user_isotope%ind > 0) deallocate(user_isotope%deuteron)

    !  n,tritium data
    if(user_isotope%int > 0) deallocate(user_isotope%tritium)

    !  transport data
    if(user_isotope%ltrn > 0) deallocate(user_isotope%transport)

    !  total data
    if(user_isotope%ltot > 0) deallocate(user_isotope%total)

    !  strpd data
    if(user_isotope%istrpd > 0) deallocate(user_isotope%strpd)

    !  scattering data
    if(user_isotope%maxterms > 0) deallocate(user_isotope%scattering)

    end subroutine

    !---------------------------------------------------------------------------------------------------
    subroutine isotxs_isotope_computeband(user_isotope)
    !  used only when need to write ISOTXS file
    type (type_isotxs_isotope), intent(INOUT) :: user_isotope
    !  local
    integer group,gprime,start_band,stop_band
    integer startlegendre,stoplegendre,ncs,iii,jjj,leg
    real(kind=sp) sum

    if(user_isotope%defined) then
        !     we need to do this for ncsmax which must be defined as specified by the isotxs data file format (i.e. idsct)
        do ncs = 1,user_isotope%nscmax
            !        figure out what it is
            iii = user_isotope%idsct(ncs)
            jjj = iii/100       ! 000,100,202,300 -> 0,1,2,3  ! the scattering type
            iii = iii - jjj*100 !   0,  0,  2,  0             ! the starting legendre order in the ncs block
            !        determine its storage position
            if(jjj == 0) startlegendre = user_isotope%total_start + iii
            if(jjj == 1) startlegendre = user_isotope%elast_start + iii
            if(jjj == 2) startlegendre = user_isotope%inela_start + iii
            if(jjj == 3) startlegendre = user_isotope%n2n_start   + iii
            stoplegendre = startlegendre                      ! lord(:)=1

            !        if we assume scattering(g,g') where g'-> g is the rule:
            !          jband becomes the length of each row in scattering
            !          ijj   is the offset from the right most non-zero column in scattering to the within group term
            !        1 2 3 4 5 6 7 8 9
            !        |---------------|     group   start_band  stop_band   jband  ijj   ds  ups
            !            |---w---|     ->    5          3          7     ->  5     3    2    2
            !        |-w---|           ->    2          1          4     ->  4     3    1    2
            !            |-------w-|   ->    7          3          8     ->  6     2    4    1
            !        the number of up   scatter groups is given as ijj-1
            !        the number of down scatter groups is given as jband-ijj

            do group = 1,user_isotope%ngroup
                start_band = group
                stop_band  = group
                searchabove: do gprime = 1,group
                    sum = 0.
                    do leg = startlegendre,stoplegendre
                        !                                                             g<-g'
                        sum = sum + abs(user_isotope%scattering(group,gprime,leg))
                    enddo
                    if(abs(sum) > 1.e-20) then
                        start_band = gprime
                        exit searchabove
                    endif
                enddo searchabove

                searchbelow: do gprime = user_isotope%ngroup,group,-1
                    sum = 0.
                    do leg = startlegendre,stoplegendre
                        !                                                             g<-g'
                        sum = sum + abs(user_isotope%scattering(group,gprime,leg))
                    enddo
                    if(abs(sum) > 1.e-20) then
                        stop_band = gprime
                        exit searchbelow
                    endif
                enddo searchbelow
                !           check the bounds on the banding
                !           if(start_band <= 0)     start_band = 1
                !           if(start_band > group)  start_band = group
                !           if(stop_band  < group)  stop_band  = group
                !           if(stop_band  > user_isotope%ngroup) stop_band  = user_isotope%ngroup
                !           set jband and ijj for this row of the scattering matrix
                user_isotope%jband(group,ncs) = stop_band - start_band + 1
                user_isotope%ijj(group,ncs)   = stop_band - group + 1
            enddo ! group
        enddo ! ncs
    endif

    end subroutine

    !---------------------------------------------------------------------------------------------------------------------------------
    subroutine isotxs_calc_group_vel(ngroup, energy, vel)
    implicit none
    integer :: ngroup
    real(kind=sp) :: energy(ngroup+1), vel(ngroup)
    ! local
    integer :: ig
    real(kind=sp) :: aneut, cvel, temp, temp1
    data aneut /9.395527e+08/
    data cvel  /2.997925e+10/
    !  V = C * sqrt(E*(2*M+E))/(M+E)
    do ig=1,ngroup
        temp =0.5*(energy(ig)+energy(ig+1))
        temp1=temp*(temp+2.*aneut)
        temp1=sqrt(temp1)
        vel(ig)=(cvel*temp1)/(aneut+temp)
    enddo
    end subroutine

    !---------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE isotxs_isotope_computeMAXORD(User_Isotope)
    IMPLICIT NONE
    ! Passed in
    TYPE (type_isotxs_isotope) :: User_Isotope

    User_Isotope%MAXORD = 0
    IF (User_Isotope%Defined) THEN
        IF (User_Isotope%LTRN          .GT. User_Isotope%MAXORD) User_Isotope%MAXORD = User_Isotope%LTRN
        IF (User_Isotope%LTOT          .GT. User_Isotope%MAXORD) User_Isotope%MAXORD = User_Isotope%LTOT
        IF (User_Isotope%TOTAL_PNTERMS .GT. User_Isotope%MAXORD) User_Isotope%MAXORD = User_Isotope%TOTAL_PNTERMS
        IF (User_Isotope%ELAST_PNTERMS .GT. User_Isotope%MAXORD) User_Isotope%MAXORD = User_Isotope%ELAST_PNTERMS
        IF (User_Isotope%INELA_PNTERMS .GT. User_Isotope%MAXORD) User_Isotope%MAXORD = User_Isotope%INELA_PNTERMS
        IF (User_Isotope%N2N_PNTERMS   .GT. User_Isotope%MAXORD) User_Isotope%MAXORD = User_Isotope%N2N_PNTERMS
        User_Isotope%MAXORD = User_Isotope%MAXORD - 1 ! Convert it from #terms to order
        IF (User_Isotope%MAXORD .LT. 0) User_Isotope%MAXORD = 0 ! Bounds check (void isotopes can have no data)
    END IF

    END SUBROUTINE

    !---------------------------------------------------------------------------------------------------------------------------------
    ! to compute maximum number of up- and down-scattering groups, for which JBAND and IJJ arrays should be defined in advance
    Subroutine isotxs_isotope_computeMAXUP(User_Isotope, MAXUP, MAXDN)
    implicit none
    type (type_isotxs_isotope) :: User_Isotope
    integer :: MAXUP, MAXDN
    integer :: isct, gg, DNG, UPG

    MAXUP = 0
    MAXDN = 0
    do isct = 1, User_Isotope%nscmax
        do gg = 1, User_Isotope%ngroup
            DNG = User_Isotope%jband(gg,isct) - User_Isotope%ijj(gg,isct)    ! number of down-scattering groups
            UPG = User_Isotope%ijj(gg,isct)-1    ! number of up-scattering groups
            MAXDN = max(MAXDN, DNG)
            MAXUP = max(MAXUP, UPG)
        end do
    end do

    End Subroutine

    !---------------------------------------------------------------------------------------------------------------------------------
    function isotxs_isotope_countrecords(isotope) result(num_records)
    implicit none
    type(type_isotxs_isotope) :: isotope
    integer :: num_records, isct
    num_records = 2
    if(isotope%ichi > 1) num_records = num_records + 1
    ! num_records = num_records + sum(2*isotope%lord(:) / (isotope%lord(:)+1))
    do isct = 1, isotope%nscmax
        ! if(isotope%lord(isct) > 0) num_records = num_records + 1
        num_records = num_records + 2*isotope%lord(isct) / (isotope%lord(isct)+1)
    end do
    end function

    !---------------------------------------------------------------------------------------------------
    Subroutine isotxs_xs_scale_uniform(isotope, weight, newiso)
    implicit none
    type(type_isotxs_isotope) :: isotope, newiso
    real(8) :: weight
    ! local
    integer :: gg, isct, jl, ju, idsct, nn, nl, jj, g_prime

    ! tran. (scalar flux weighted), total, ngamma
    newiso%transport(:,:) = isotope%transport(:,:) * weight
    newiso%total(:,:) = isotope%total(:,:) * weight
    newiso%ngamma(:) = isotope%ngamma(:) * weight
    ! fission, nu
    if(newiso%ifis > 0) newiso%fission(:) = isotope%fission(:) * weight

    ! chi - no change

    ! others
    if(newiso%ialf > 0) newiso%alpha(:)    = isotope%alpha(:) * weight
    if(newiso%inp > 0)  newiso%proton(:)   = isotope%proton(:) * weight
    if(newiso%in2n > 0) newiso%n2n(:)      = isotope%n2n(:) * weight
    if(newiso%ind > 0)  newiso%deuteron(:) = isotope%deuteron(:) * weight
    if(newiso%int > 0)  newiso%tritium(:)  = isotope%tritium(:) * weight
    if(newiso%istrpd > 0) then    ! don't think it is really used any more
        do jj = 1, newiso%istrpd
            newiso%strpd(:,jj) = isotope%strpd(:,jj) * weight    ! should be weighted by current
        end do
    end if

    ! scattering
    do isct = 1, isotope%nscmax
        idsct = isotope%idsct(isct) / 100
        nn = isotope%idsct(isct) - idsct*100
        do gg = 1, isotope%ngroup
            if(isotope%jband(gg,isct) == 0) cycle
            ! determine range of groups that scatter into this group (gg)
            jl = gg + isotope%ijj(gg,isct) - isotope%jband(gg,isct)
            ju = gg + isotope%ijj(gg,isct) - 1
            do nl = 1, isotope%lord(isct)
                ! determine scattering block storage location (note we don't support subblock)
                select case (idsct)
                case (0)
                    jj = isotope%total_start + nn + nl - 1
                case (1)
                    jj = isotope%elast_start + nn + nl - 1
                case (2)
                    jj = isotope%inela_start + nn + nl - 1
                case (3)
                    jj = isotope%n2n_start + nn + nl - 1
                end select
                ! scale
                do g_prime = jl, ju
                    newiso%scattering(gg,g_prime,jj) = isotope%scattering(gg,g_prime,jj) * weight
                end do
            end do
        end do
    end do  ! loop over scattering blocks

    End Subroutine

    !---------------------------------------------------------------------------------------------------
    Subroutine isotxs_xs_scale_groupwise(isotope, ngroup, lord, weight, newiso)
    implicit none
    type(type_isotxs_isotope) :: isotope, newiso
    integer :: ngroup, lord
    real(8) :: weight(ngroup, lord)
    ! local
    integer :: isct, jl, ju, jj, gg, g_prime, nn, nl, ilord, idsct

    if(ngroup /= isotope%ngroup) then
        write(0,'("[isotxs_xs_scale_groupwise]...Error: inconsistent energy groups of weighting function to scale XS of isotopes")')
        call abort('')
    end if

    ! tran. (scalar flux weighted), total, ngamma
    newiso%transport(:,1) = isotope%transport(:,1) * weight(:,1)
    newiso%total(:,1) = isotope%total(:,1) * weight(:,1)
    newiso%ngamma(:) = isotope%ngamma(:) * weight(:,1)
    ! fission, nu
    if(newiso%ifis > 0) newiso%fission(:) = isotope%fission(:) * weight(:,1)

    ! chi - no change

    ! others
    if(newiso%ialf > 0) newiso%alpha(:)    = isotope%alpha(:) * weight(:,1)
    if(newiso%inp > 0)  newiso%proton(:)   = isotope%proton(:) * weight(:,1)
    if(newiso%in2n > 0) newiso%n2n(:)      = isotope%n2n(:) * weight(:,1)
    if(newiso%ind > 0)  newiso%deuteron(:) = isotope%deuteron(:) * weight(:,1)
    if(newiso%int > 0)  newiso%tritium(:)  = isotope%tritium(:) * weight(:,1)

    ! scattering
    do isct = 1, isotope%nscmax
        idsct = isotope%idsct(isct) / 100
        nn = isotope%idsct(isct) - idsct*100
        ilord = nn + 1
        !
        ! use scalar flux weight if no higher-order flux moments, check this later
        !
        if(ilord > lord) ilord = 1
        do gg = 1, ngroup
            if(isotope%jband(gg,isct) == 0) cycle
            ! determine range of groups that scatter into this group (gg)
            jl = gg + isotope%ijj(gg,isct) - isotope%jband(gg,isct)
            ju = gg + isotope%ijj(gg,isct) - 1
            do nl = 1, isotope%lord(isct)
                ! determine scattering block storage location (note we don't support subblock)
                select case (idsct)
                case (0)
                    jj = isotope%total_start + nn + nl - 1
                case (1)
                    jj = isotope%elast_start + nn + nl - 1
                case (2)
                    jj = isotope%inela_start + nn + nl - 1
                case (3)
                    jj = isotope%n2n_start + nn + nl - 1
                end select
                ! scale
                do g_prime = jl, ju
                    newiso%scattering(gg,g_prime,jj) = isotope%scattering(gg,g_prime,jj) * weight(g_prime,ilord)
                end do
            end do
        end do
    end do  ! loop over scattering blocks

    ! high-order total and transport xs
    if(newiso%istrpd > 0) then    ! don't think it is really used any more
        do jj = 1, newiso%istrpd
            newiso%strpd(:,jj) = isotope%strpd(:,jj) * weight(:,1)    ! should be weighted by current
        end do
    end if

    ! -- incomplete ...

    End Subroutine

    !---------------------------------------------------------------------------------------------------
    ! compute total scattering production matrices and store total matrices only 
    !Subroutine isotxs_convert_scatmx_to_total(isotxs, isotope)
    !    implicit none
    !    type(type_isotxs) :: isotxs
    !    type(type_isotxs_isotope) :: isotope
    !End Subroutine
        


    !---------------------------------------------------------------------------------------------------
    ! get total scattering (production) matrices, multiplicity of (n,2n) scattering is included
    Subroutine isotxs_get_total_scatmx(isotope, ngroup, maxord, scat)
    implicit none
    type(type_isotxs_isotope) :: isotope
    integer :: ngroup, maxord
    real(kind=sp), intent(out) :: scat(ngroup,ngroup,0:maxord)
    ! local
    integer :: gg, l, n
    real(kind=sp) :: n2nfactor
    logical :: total_scattering_available = .false.

100 format('[isotxs_get_total_scatmx]...Error: ',A,1X,I4)
200 format('[isotxs_get_total_scatmx]...Warning: ',A,1X,I4)

    if(isotope%total_pnterms>0) then  ! total scattering matrices exist
        if(isotope%total_pnterms >= maxord+1) then
            n = isotope%total_start
            do l = 0, maxord
                scat(:,:,l) = isotope%scattering(:,:,n+l)
            end do
            total_scattering_available = .true.
        else
            if(isotope%elast_start<=0 .or. isotope%elast_pnterms<maxord+1) then
                write(0,100) 'requested scattering matrices of '//isotope%hisonm//' contain higher Legendre orders than available in isotxs'
                call abort('')
            end if
        end if
    end if

    if(total_scattering_available) return

    n = isotope%elast_start
    do l = 0, maxord
        scat(:,:,l) = isotope%scattering(:,:,n+l)
    end do

    if(isotope%inela_pnterms > 0) then
        n = isotope%inela_start
        do l = 0, isotope%inela_pnterms-1
            scat(:,:,l) = scat(:,:,l) + isotope%scattering(:,:,n+l)
        end do
    else 
        write(0,200) 'while total scattering matrix not available, inelastic scattering matrix of '//isotope%hisonm//' is not (explicitly) included in isotxs'
    end if

    if (isotope%n2n_pnterms>0) then
        n2nfactor = isotxs_Isotope_n2nfactor(isotope)   ! = 1.0 / 2.0  for production / reaction based n2n matrix
        n = isotope%n2n_start
        do l = 0, isotope%n2n_pnterms-1
            scat(:,:,l) = scat(:,:,l) + n2nfactor * isotope%scattering(:,:,n+l)
        end do
    else
        write(0,200) 'while total scattering matrix not available, n2n scattering matrix of '//isotope%hisonm//' is not (explicitly) included in isotxs'
    end if

    End Subroutine

    !---------------------------------------------------------------------------------------------------
    ! get n2n factor (= 1.0 / 2.0 for production / reaction based n2n scattering matrix)
    real(kind=sp) Function isotxs_isotope_n2nfactor(isotope)
    implicit none
    type (type_isotxs_isotope) :: isotope
    ! local
    integer :: Itemp,Group,Gprime
    real(kind=sp) :: Sum_of_principle_n2n,Sum_of_scattering_n2n,Return_N2N

    Return_N2N = 0.0d0 ! Return zero if nothing can be said
    IF (isotope%defined) THEN
        IF ((isotope%N2N_PNTERMS .NE. 0) .AND. (isotope%IN2N .NE. 0)) THEN
            Sum_of_principle_n2n  = 0.0d0
            Sum_of_scattering_n2n = 0.0d0
            Itemp = isotope%N2N_START
            DO Group = 1,isotope%NGROUP
                Sum_of_principle_n2n = Sum_of_principle_n2n + isotope%N2N(Group)
                DO Gprime = 1,isotope%NGROUP
                    Sum_of_scattering_n2n = Sum_of_scattering_n2n + isotope%SCATTERING(Gprime,Group,Itemp)
                END DO
            END DO
            IF (Sum_of_principle_n2n / Sum_of_scattering_n2n .GT. 0.75D0) THEN
                Return_N2N = 2.0D0 ! n2n scattering matrix is reaction based
            ELSE
                Return_N2N = 1.0D0 ! n2n scattering matrix is production based
            END IF
        ELSE IF (isotope%IN2N .NE. 0) THEN
            Return_N2N = 1.0D0 ! It is safe to assume that the user used production based n2n scattering matrices and n2n reaction principle as ISOTXS specifies
        END IF
    END IF

    isotxs_Isotope_n2nfactor = Return_N2N
    
    end function

    !---------------------------------------------------------------------------------------------------
    !subroutine isotxs_check_xs_balance(isotope)
    !type(type_isotxs_isotope) :: isotope
    !
    !end subroutine
    
    !---------------------------------------------------------------------------------------------------
    ! using high-order total XSs to make transport correction to self-scattering XSs
    ! only applied to ISOTXS that contains high-order total cross sections
    subroutine isotxs_transport_scatgg(isotope)
    type(type_isotxs_isotope) :: isotope
    ! local
    integer :: isct, idsct, nn, nl, jj, gg

    100 format('[isotxs_transport_scatgg]...',A,1X,I4)
    if(isotope%total_pnterms > isotope%ltot .or. isotope%elast_pnterms > isotope%ltot) then
        write(0,100) 'Error: missing some high-order total cross sections'
        call abort('')
    end if
    if(isotope%scatgg_corrected) then
        write(0,100) 'Error: extended transport approximation has been made to this isotope '//isotope%hisonm
        call abort('')
    end if

    ! Sigs(L,g->g) <= Sigs(L,g->g) + (Sigt(0,g) - Sigt(L,g))
    do isct = 1, isotope%nscmax
        idsct = isotope%idsct(isct) / 100
        if(idsct > 1) cycle     ! correct total/elastic scattering only, assumed that total scattering and component scattering matrices not occur simultaneously
        nn = isotope%idsct(isct) - idsct*100
        if(nn == 0) cycle
        do nl = 1, isotope%lord(isct)
            if(idsct == 0) then
                jj = isotope%total_start + nn + nl - 1
            else ! idsct == 1
                jj = isotope%elast_start + nn + nl - 1
            end if
            do gg = 1, isotope%ngroup
                isotope%scattering(gg,gg,jj) = isotope%scattering(gg,gg,jj) + isotope%total(gg,1) - isotope%total(gg,nn+1)
            end do
        end do
    end do
    isotope%scatgg_corrected = .true.

    end subroutine

    !---------------------------------------------------------------------------------------------------
    ! using high-order total XSs to make transport correction to self-scattering XSs
    ! only applied to ISOTXS that contains high-order total cross sections
    subroutine isotxs_revert_transport_scatgg(isotope)
    type(type_isotxs_isotope) :: isotope
    ! local
    integer :: isct, idsct, nn, nl, jj, gg
    
    100 format('[isotxs_revert_transport_scatgg]...',A,1X,I4)
    if(isotope%total_pnterms > isotope%ltot .or. isotope%elast_pnterms > isotope%ltot) then
        write(0,100) 'Error: missing some high-order total cross sections'
        call abort('')
    end if
    if(.not. isotope%scatgg_corrected) return
    
    ! Sigs(L,g->g) <= Sigs(L,g->g) - (Sigt(0,g) - Sigt(L,g))
    do isct = 1, isotope%nscmax
        idsct = isotope%idsct(isct) / 100
        if(idsct > 1) cycle     ! correct total/elastic scattering only, assumed that total scattering and component scattering matrices not occur simultaneously
        nn = isotope%idsct(isct) - idsct*100
        if(nn == 0) cycle
        do nl = 1, isotope%lord(isct)
            if(idsct == 0) then
                jj = isotope%total_start + nn + nl - 1
            else ! idsct == 1
                jj = isotope%elast_start + nn + nl - 1
            end if
            do gg = 1, isotope%ngroup
                isotope%scattering(gg,gg,jj) = isotope%scattering(gg,gg,jj) - isotope%total(gg,1) + isotope%total(gg,nn+1)
            end do
        end do
    end do
    isotope%scatgg_corrected = .false.

    end subroutine

    end module
