Subroutine isotxs_read_file(lib_file, isotxs, NOUT, iread_option)
#include 'Custom_Types.h'    
    use m_isotxs
    use m_kind
    use NE_Kind, only: NE_Kind_GetFreeLogicalUnit
    implicit none 
    character(*) :: lib_file
    integer :: NOUT
    integer :: iread_option   ! 0 / 1  read header records / all data
    type (type_isotxs) :: isotxs
    ! local
    integer :: NIN, IOS, ii, jj, kk
    integer :: niso, ngroup, nscmax, isize_scratch, isize_scratch2
    type (type_isotxs_isotope), pointer :: isotope
    real(kind=sp), allocatable :: scratch(:,:), scratch2(:)
    
    
100 format('[isotxs_read_file]...Error: ',A)

    NIN = NE_Kind_GetFreeLogicalUnit()
    open(NIN,file=lib_file,status='old',form='unformatted',iostat=IOS)
    if(IOS .NE. 0) then 
        write(NOUT,100) 'cannot open cross section library <'//trim(lib_file)//'>'
        call abort('') 
    end if
    isotxs%file_unit = NIN
    
    ! read first 4 records (0D ~ 3D records)
    call isotxs_read_header(NIN, isotxs, 1)
    isotxs%nskip_header = 3
    if (isotxs%ichist > 1) isotxs%nskip_header = 4
    if (iread_option == 0) then
        rewind(NIN)
        return
	end if
    
    niso = isotxs%niso
    ngroup = isotxs%ngroup
    nscmax = isotxs%nscmax

    isize_scratch = 7 + 2*(isotxs%maxord+1) + 6       ! for principal XS
    ! for scattering XS, (maxord+1) elastic + 2 inelastic + 1 n2n matrices at most, 
    !   it is assumed that total scattering matrix would not occur with component matrices at same time
    isize_scratch2 = ngroup * ngroup * (isotxs%maxord+4)  
    
    allocate(isotxs%isotopedata(niso), scratch(ngroup, isize_scratch), scratch2(isize_scratch2), stat=IOS)
    if(IOS .NE. 0) then 
        write(NOUT,100) 'failed to allocate storage for reading cross sections'
        call abort('') 
    end if
    
    ! read all isotopes
    do ii = 1, niso
        isotope => isotxs%isotopedata(ii)
        ! allocate for isotope header
        isotope%hisonm = isotxs%hisonm(ii)
        !write(0,*) 'reading isotope '//isotope%hisonm
        isotope%nscmax = nscmax
        isotope%ngroup = ngroup
        allocate(isotope%idsct(nscmax), isotope%lord(nscmax), isotope%jband(ngroup,nscmax), isotope%ijj(ngroup,nscmax))
        ! read 4D record and allocate space for this isotope
        call isotxs_read_isotope(NIN, ngroup, isotope, isize_scratch, isize_scratch2, scratch, scratch2, 1)
        call isotxs_allocate_isotope_xs(isotope, ngroup, nscmax, isotope%lord, isotope%idsct)
        isotope%defined = .true.
        ! read 5D through 7D records of this isotope
        call isotxs_read_isotope(NIN, ngroup, isotope, isize_scratch, isize_scratch2, scratch, scratch2, 2)
    end do
    isotxs%defined = .true.
    
End Subroutine