Subroutine reload_macxs(NOUT, ixs)
    use m_reload
    use COMPXS_IO
    use m_isotxs
    implicit none 
    integer :: NOUT, ixs  ! ixs = 1 / 2, compute only heating xs / all XSs, i.e. COMPXS
    ! local
    type (type_isotxs) :: isotxs
    type (type_isotxs_isotope), pointer :: isotope
    integer :: NGROUP, IOS, ii, jj, icmp, jcmp
    character(len=8) :: ANAME
    real(4), allocatable :: xs(:)
    real(8) :: fdum

    100 format('[reload_macxs]...',A,1X,I6)

    call isotxs_read_file(isotxs_file, isotxs, NOUT, 1)
    ngroup = isotxs%ngroup
    kerma%ncmp = num_reload_comp
    kerma%ngroup = ngroup
    ! save group structure
    allocate(grpstr(ngroup+1))
    grpstr = isotxs%energy
    
    if(gamma_heating) then
        write(NOUT,100) 'Error: not implemented for gamma heating yet'
        call abort
    end if

    if(ixs == 1) then
        allocate(kerma%nheat(ngroup,num_reload_comp+num_test_comp), kerma%gheat(0,num_reload_comp+num_test_comp))
        kerma%nheat = 0.0D0
        IOS = 0
        do icmp = 1, num_reload_comp
            do ii = 1, RFUEL(icmp)%niso
                ANAME = RFUEL(icmp)%isonam(ii)
                do jj = 1, isotxs%niso
                    if(ANAME == isotxs%hisonm(jj)) then
                        isotope => isotxs%isotopedata(jj)
                        exit
                    end if
                end do
                if (jj > isotxs%niso) then
                    write(NOUT,100) 'Error: isotope '//ANAME//' in reloaded composition '//RFUEL(icmp)%name//' is not found in ISOTXS'
                    IOS = 1
                else
                    xs = isotope%ngamma
                    if(isotope%ialf > 0) xs = xs + isotope%alpha
                    if(isotope%inp > 0) xs = xs + isotope%proton
                    if(isotope%ind > 0) xs = xs + isotope%deuteron
                    if(isotope%int > 0) xs = xs + isotope%tritium
                    xs = xs * isotope%ecapt
                    if(isotope%ifis > 0) xs = xs + isotope%efiss * isotope%fission
                    kerma%nheat(:,icmp) = kerma%nheat(:,icmp) + RFUEL(icmp)%density(ii) * xs
                end if
            end do
        end do
        if (measure_uncertainty) then
            do jcmp = 1, num_test_comp
                icmp = num_reload_comp + jcmp
                do ii = 1, TFUEL(jcmp)%niso
                    ANAME = TFUEL(jcmp)%isonam(ii)
                    do jj = 1, isotxs%niso
                        if(ANAME == isotxs%hisonm(jj)) then
                            isotope => isotxs%isotopedata(jj)
                            exit
                        end if
                    end do
                    if (jj > isotxs%niso) then
                        write(NOUT,100) 'Error: isotope '//ANAME//' in test composition '//TFUEL(jcmp)%name//' is not found in ISOTXS'
                        IOS = 1
                    else
                        xs = isotope%ngamma
                        if(isotope%ialf > 0) xs = xs + isotope%alpha
                        if(isotope%inp > 0) xs = xs + isotope%proton
                        if(isotope%ind > 0) xs = xs + isotope%deuteron
                        if(isotope%int > 0) xs = xs + isotope%tritium
                        xs = xs * isotope%ecapt
                        if(isotope%ifis > 0) xs = xs + isotope%efiss * isotope%fission
                        kerma%nheat(:,icmp) = kerma%nheat(:,icmp) + TFUEL(jcmp)%density(ii) * xs
                    end if
                end do
            end do
        end if
        if (IOS > 0) call abort
        kerma%defined = .true.
        do ii = 1, isotxs%niso
            call isotxs_deallocate_isotope(isotxs%isotopedata(ii))
        end do
        call isotxs_deallocate_header(isotxs)
    else 
        write(NOUT,100) 'Error: generation of COMPXS for reloaded fuel is disabled until internal perturbation calculation is implemented'
        call abort
        allocate(kerma%nheat(ngroup,num_reload_comp), kerma%gheat(0,num_reload_comp), xs(ngroup))
        ! we disable generation of COMPXS for reloaded composition unless we do perturbation calculation internally (w/o PERSENT)

        !call COMPXS_ASSIGNPRINTINFO(NOUT)
        !! load in existing COMPXS
        !call COMPXS_IMPORT(NOUT, COMPXS, 'COMPXS')
        !! compute macro XS for reloaded fuel
        !call COMPXS_DEFINE(RFUELXS, num_reload_comp, ngroup, 0, 0, isotxs%maxord)
        !call COMPXS_ZERO(RFUELXS)
        !do icmp = 1, num_reload_comp
        !    call isotxs_add_zone_to_compxs(NOUT, RFUELXS, icmp, RFUEL(icmp)%niso, &
        !        RFUEL(icmp)%isonam, RFUEL(icmp)%density, isotxs)
        !end do
        !kerma%nheat = RFUELXS%PC
        !kerma%defined = .true.
    end if

End Subroutine
