Subroutine reload_input(NOUT)
#include 'Custom_Types.h'
    use NE_Kind
    use m_reload
    implicit none
    integer :: NOUT
    logical :: yes
    ! &control
    namelist /control/ flux_data, isotxs_file, gamma_heating, consider_keff, consider_worth, use_persent, consider_burnup, &
        adjust_cycle_length, num_reload_comp, fuel_regions, num_fuel_region, num_hex_sector, burnup_data, &
        fixed_peak_power, peak_power_limit, assembly_power_limit, reload_refpwr, update_refpwr, &
        objective_func, min_burnup, min_deltak, max_deltak, max_worth, dif3d_check, dif3d_check_all, normalize_power, &
        min_cycle_length, target_eoc_keff, deprate_ratio, eps_cycle, &
        null_refueling, num_persent_case, candidate_range, candidate_range_size, measure_uncertainty, num_test_comp, fixed_sequence, &
        margin_power_limit, eps_power, wgt_pekfac
    ! &material 
    type isotope_info_type
        Char8 :: name = 'none'   ! unique name
        Real4 :: density  = 0.0E0    ! atom density
    end type
    Char8, allocatable :: composition_name(:)
    type (isotope_info_type), allocatable :: composition(:,:)
    Char8, allocatable :: testcomp_name(:)
    type (isotope_info_type), allocatable :: testcomp(:,:)
    namelist /material/ composition_name, composition, testcomp_name, testcomp
    ! &
    ! local
    integer :: NIN, IOS, ii, jj, NISO, NCMP, kk = 0

    100 format('[reload_input]...',A,1X,I6)

    NIN = NE_Kind_GetFreeLogicalUnit()
    open(NIN, file='reload.in', status='old', form='formatted', iostat=IOS)
    if(IOS /= 0) then
        write(NOUT,100) 'Error: cannot open input file "reload.in"'
        call abort
    end if

    ! control options
    read(NIN,nml=control, iostat=IOS)
    if(IOS /= 0) then
        write(NOUT,100) 'Error: failed to read input block &control'
        call abort
    end if
    call NE_Kind_UpperCaseString(flux_data)
    IOS = 0
    if(gamma_heating) then
        write(NOUT,100) 'Error: gamma heating is not considered yet'
        IOS = 1
    end if
    if(num_hex_sector < 0 .or. (num_hex_sector > 2 .and. num_hex_sector /=6)) then
        write(NOUT,100) 'Error: valid num_hex_sector value = 0, 1, 2, or 6'
        IOS = 1
    end if
    call read_cycle_length(NOUT, 'ABURN', cycle_length)
    if(adjust_cycle_length) then
        if(.not. consider_keff) then
            consider_keff = .true.
            min_deltak = 0.000001
            write(NOUT,100) 'Note: reactivity change is considered when cycle length is adjustable and min_deltak will &
            &be determined by min_cycle_length'
        end if
        consider_burnup = .true.    ! always consider discharge burnup when cycle length is adjustable
        if(min_cycle_length < 0.0) then
            write(NOUT,100) 'Error: invalid specification of minimum cycle length, should be positive'
            IOS = 1
        end if
        if(cycle_length / min_cycle_length - 1.0d0 < -eps_cycle) &
        &    write(NOUT,100) 'Warning: cycle length of this fuel cycle is smaller than the specified minimum cycle length'
        if(target_eoc_keff < 0.0) then
            write(NOUT,100) 'Error: invalid specification of targeted keff at end of cycle, should be positive and typically ~ 1.0'
            IOS = 1
        end if
    end if
    do ii = 1, num_objectives
        if(objective_func == defined_objectives(ii)) exit
    end do
    if(ii > num_objectives) then 
        write(NOUT,100) 'Error: invalid option for objective_func -> '//objective_func
        write(NOUT,'(7X,"valid options are: ",(5(A16,1X)))') defined_objectives
        IOS = 1
    end if
    if(objective_func == 'burnup') then
        consider_burnup = .true.
    else if(objective_func == 'reactivity') then
        use_persent = .true.
    else if(objective_func == 'reactivity_slow') then
        use_persent = .true.
    else if(objective_func == 'power_shape' .or. objective_func == 'power_deep') then
        if(reload_refpwr == 'flat') update_refpwr = .false.
    end if
    if(consider_worth) then
        use_persent = .true.
        if(max_worth < 0.0) max_worth = 0.002D0   ! 200 pcm
    end if
    if(consider_keff) then
        if(.not. use_persent) dif3d_check = .true.
        if(max_deltak < 0.0) then
            if(consider_worth) then
                max_deltak = max_worth
            else
                max_deltak = 0.002D0   ! 200 pcm
            end if
        else
            if(consider_worth .and. max_deltak > max_worth) then
                write(NOUT,100) 'Error: invalid specification of reactivity constraints, max_deltak > max_worth'
                IOS = 1
            end if
        end if
        if(min_deltak > max_deltak) then
            write(NOUT,100) 'Error: invalid specification of reactivity constraints, min_deltak > max_deltak'
            IOS = 1
        end if
    end if
    if(consider_burnup) then
        inquire(file=burnup_data, exist=yes)
        if(.not. yes) then
            write(NOUT,100) 'Error: burnup dataset <'//trim(burnup_data)//'> does not exist'
            IOS = 1
        end if
    end if
    if(null_refueling) then
        adjust_cycle_length = .false.
        consider_worth = .true.
        use_persent = .true.
    end if
    if(num_fuel_region > max_num_fuel_regions) then
        write(NOUT,100) 'Error: number of fuel regions exceeds the allowed max_num_fuel_regions,', max_num_fuel_regions
        write(NOUT,100) '       please modify m_reload.f90'
        IOS = 1
    end if
    if(IOS > 0) call abort

    ! reloaded fuel compositions
    NCMP = num_reload_comp
    if(consider_worth) NCMP = num_reload_comp + 1
    allocate(composition_name(NCMP), composition(max_niso_per_comp, NCMP), RFUEL(NCMP),   stat=IOS)
    IF (measure_uncertainty) allocate(testcomp_name(num_test_comp), testcomp(max_niso_per_comp, num_test_comp), TFUEL(num_test_comp),   stat=IOS)
    if(IOS /= 0) then
        write(NOUT,100) 'Error: failed to allocate storage for reloaded compositions'
        call abort
    end if
    !  read material information -----------------------------------------------------------------------------------
    read(NIN, nml=material, iostat=IOS)
    if(IOS /= 0) then
        write(NOUT,100) 'Error: failed to read input block &material',IOS
        call abort
    end if

    do jj = 1, NCMP
        do ii = 1, max_niso_per_comp
            if(composition(ii,jj)%name == 'none') then
                NISO = ii - 1
                exit
            end if
        end do
        if(ii == 1) then
            write(NOUT,100) 'Error: no isotopes in composition '//composition_name(jj)
            call abort 
        else if(ii > max_niso_per_comp) then
            write(NOUT,100) 'Error: number of isotopes in '//composition_name(jj)//' exceeds ',max_niso_per_comp
            write(NOUT,100) '       increase max_niso_per_comp in m_reload.f90'
            call abort
        end if
        RFUEL(jj)%name = composition_name(jj)
        RFUEL(jj)%niso = NISO
        allocate(RFUEL(jj)%isonam(NISO), RFUEL(jj)%density(NISO))
        do ii = 1, NISO
            RFUEL(jj)%isonam(ii) = composition(ii,jj)%name
            RFUEL(jj)%density(ii) = composition(ii,jj)%density
        end do
        RFUEL(jj)%defined = .true.
    end do
    ! variables for uncertainty
    if (measure_uncertainty) then
        NCMP = num_test_comp
        do jj = 1, NCMP
            do ii = 1, max_niso_per_comp
                if(testcomp(ii,jj)%name == 'none') then
                    NISO = ii - 1
                    exit
                end if
            end do
            if(ii == 1) then
                write(NOUT,100) 'Error: no isotopes in testcomp '//testcomp_name(jj)
                call abort 
            else if(ii > max_niso_per_comp) then
                write(NOUT,100) 'Error: number of isotopes in '//testcomp_name(jj)//' exceeds ',max_niso_per_comp
                write(NOUT,100) '       increase max_niso_per_comp in m_reload.f90'
                call abort
            end if
            TFUEL(jj)%name = testcomp_name(jj)
            TFUEL(jj)%niso = NISO
            allocate(TFUEL(jj)%isonam(NISO), TFUEL(jj)%density(NISO))
            do ii = 1, NISO
                TFUEL(jj)%isonam(ii) = testcomp(ii,jj)%name
                TFUEL(jj)%density(ii) = testcomp(ii,jj)%density
            end do
            TFUEL(jj)%defined = .true.
        end do
    end if
    close(NIN)

    ! positions of previously reloaded assemblies (just a few recent cycles) for distance check 
    inquire(file=reload_log, exist=yes)
    if(yes) then
        open(NIN, file=reload_log, status='old', form='formatted')
        read(NIN,*) num_previous_site, disperse_mindist
        disperse_reload = (num_previous_site > 0)
        if(disperse_reload) then
            allocate(previous_pos(num_previous_site))
            read(NIN,*) previous_pos(1:num_previous_site)
        end if
        close(NIN)
    else
        disperse_reload = .false.
    end if
    
    ! check external file for manual setting of refueling
    inquire(file=reload_manual, exist=yes)
    if(yes) then
        open(NIN, file=reload_manual, status='old', form='formatted')
        read(NIN,*, iostat=IOS) manual_reload_pos(1:2)
        if(IOS /= 0) then
            write(NOUT,100) reload_manual//' Error: cannot read the ring and position number for manually set refueling position'
            call abort
        end if
        manual_reload = .true.
        !use_persent   = .false. ! for now, disable it unless we need to calculate total worth later
        dif3d_check   = .true. 
        objective_func = 'peak_power'
    else
        manual_reload = .false.
    end if

End Subroutine
    
! ----------------------------------------------------------------------
! read cycle length from interface file ABURN
Subroutine read_cycle_length(NOUT, filename, cycle_length)
    use NE_Kind, only: NE_Kind_GetFreeLogicalUnit
    use m_utility
    implicit none 
    integer :: NOUT
    character(len=*) :: filename
    real(8) :: cycle_length
    ! local
    integer :: NIN, nskip, MCARD, ncards(99), freeform
    real(8) :: fdum
    character(8) :: ANAME
    character(72) :: line

    100 format('[read_cycle_length]...',A,1X,I8)

    NIN = NE_Kind_GetFreeLogicalUnit()
    open(NIN, file=filename, status='old', form='formatted')
    call scan_arc_file_header(NIN,MCARD,ncards,freeform,ANAME)
    if(ANAME /= 'A.BURN  ') then 
        write(NOUT,100) 'Error: the opened interface file is not A.BURN'
        call abort
    end if
    if(MCARD < 3) then
        write(NOUT,100) 'Error: this A.BURN file does not contain cycle length information (type 03 card)'
        call abort
    end if

    nskip = sum(ncards(1:2))
    call SkipLines(NIN,nskip)
    if(freeform == 1) then
        read(NIN,*) nskip, nskip, fdum, fdum, cycle_length
    else
        read(NIN,'(A)') line 
        read(line(37:48),*) cycle_length
    end if
    close(NIN)

End Subroutine