Subroutine reload_power(NOUT)
    use m_reload
    use NE_Kind, only: NE_Kind_GetFreeLogicalUnit
    implicit none
    integer :: NOUT
    ! local
    integer :: NINTI, NINTJ, NINTK, UNIT
    integer :: IASS, IZ, IX, IY, ii, jj, iload
    real(8) :: AREA, normfac, FDUM
    real(8), parameter :: sqrt3 = 1.732050807568877D0
    logical :: refpwr_exist, write_refpwr
    integer :: peak_ass, peak_pos(2), peak_ass2, peak_pos2(2)
    integer :: local_peak_nod, new_peak_ass, new_peak_pos(2)
    real(8) :: peak_pwr, peak_pwr2, peak_den, peak_den2 ! peak assembly power (W), peak power density (W/cc)
    real(8) :: local_peak_den, new_peak_pwr, new_peak_den

    100 format('[reload_power]...',A,1X,I6)
    NINTI = GEODST%NINTI
    NINTJ = GEODST%NINTJ
    NINTK = GEODST%NINTK

    ! compute node volumes
    allocate(node_vol(NINTK))
    if(GEODST%IGOM==10 .or. GEODST%IGOM==18) then   ! hexagonal
        AREA = GEODST%XMESH(2) - GEODST%XMESH(1)
        AREA = sqrt3*0.5*AREA*AREA
    else if(GEODST%IGOM==6 .or. GEODST%IGOM==14) then ! cartesian
        AREA = (GEODST%XMESH(2) - GEODST%XMESH(1)) * (GEODST%YMESH(2) - GEODST%YMESH(1))
    else
        write(NOUT,100) 'Error: only hexagonal/Cartesian geometries are supported'
        call abort
    end if
    if(GEODST%IGOM==14 .or. GEODST%IGOM==18) then
        do ii = 1, GEODST%NINTK
            node_vol(ii) = AREA * (GEODST%ZMESH(ii+1) - GEODST%ZMESH(ii))
        end do
    else    ! 2D cases
        node_vol(1) = AREA
    end if

    ! read fuel assembly powers at BOC
    if(allocated(bocpwr)) deallocate(bocpwr)
    if(allocated(bocden)) deallocate(bocden)
    allocate(bocpwr(num_xynod), bocden(num_xynod), peknod(num_xynod))
    call get_peak_power(boc_PWDINT, NOUT, bocpwr, peak_ass, peak_pwr, bocden, peknod, peak_pos, peak_den)
    boc_peak_assembly_power = peak_pwr
    boc_peak_power_density = peak_den
    write(NOUT,'(/,20X,"*** Assembly powers at the beginning of cycle (BOC) ***")')
    if(num_hex_sector > 0) then
        write(NOUT,'(" IRING = hexagon ring number, IPOS = position number, Node index starts from bottom to top",/, &
        &"Index  IRING  IPOS   Power (W)    Peak density (W/cc)    Peak node",/)')
    else
        write(NOUT,'(" IX = assembly position along X dimension, IY = assembly position along Y dimension",/, &
        &"Index    IX    IY    Power (W)    Peak density (W/cc)    Peak node",/)')
    end if
    do ii = 1, num_potential_site
        IASS = fuel_pos(ii)
        if(num_hex_sector > 0) then
            write(NOUT,'(I5,2I6,2ES15.5,7X,I6)') IASS, hexagon_rp(1:2,IASS), bocpwr(IASS), bocden(IASS), peknod(IASS)
        else
            call reload_ass2xy(IASS, num_xnode, IX, IY)
            write(NOUT,'(I5,2I6,2ES15.5,7X,I6)') IASS, IX, IY, bocpwr(IASS), bocden(IASS), peknod(IASS)
        end if
    end do
    write(NOUT,'(/,"Peak assembly power at BOC occurs in assembly ",I6,", power (watts) =",ES13.5)') peak_ass, peak_pwr
    write(NOUT,'("Peak node-averaged power density at BOC occurs in assembly ",I6,", axial node ",I6,", density (w/cc) =",ES13.5,/)') &
    &   peak_pos(1:2), peak_den

    ! get reference assembly powers
    allocate(refpwr(num_xynod))
    refpwr(:) = 0.0d0
    write_refpwr = .false.
    if(reload_refpwr == 'flat') then
        fdum = PWDINT%POWER / num_potential_site    ! average assembly power
        ! it is ok to assume all power generated in fuel assemblies since the average value is used as a constant scaling factor
        do ii = 1, num_potential_site
            IASS = fuel_pos(ii)
            refpwr(IASS) = fdum
        end do
    else    ! 'bol' or 'user'
        inquire(file=ref_power, exist=refpwr_exist)
        UNIT = NE_Kind_GetFreeLogicalUnit()
        if(refpwr_exist) then
            open(unit=UNIT, file=ref_power, status='old', form='formatted')
            do ii = 1, num_potential_site
                read(UNIT,*, iostat=jj) IASS, refpwr(IASS)
                if(jj /= 0) then
                    write(NOUT,100) 'Error: I/O error encountered when reading reference assembly powers from '//ref_power
                    call abort 
                end if
            end do
            close(UNIT)
            if(update_refpwr .and. boc_peak_assembly_power < maxval(refpwr)) then
                write_refpwr = .true.
                refpwr = bocpwr
            end if
        else
            refpwr = bocpwr
            write_refpwr = .true.
        end if
        ! record reference power distribution
        if(write_refpwr) then
            open(unit=UNIT, file=ref_power, status='replace', form='formatted')
            do ii = 1, num_potential_site
                IASS = fuel_pos(ii)
                write(UNIT,'(I5,2X,ES13.5)') IASS, refpwr(IASS)                
            end do
            close(UNIT)
        end if
    end if

    ! compute peak assembly power at EOC before reloading
    if(allocated(eocpwr)) deallocate(eocpwr)
    if(allocated(eocden)) deallocate(eocden)
    if(allocated(asspwr)) deallocate(asspwr)
    if(allocated(pekpwr)) deallocate(pekpwr)
    if(allocated(pekfac)) deallocate(pekfac)
    if(allocated(pekass)) deallocate(pekass)
    if(allocated(pekden)) deallocate(pekden)
    if(allocated(pekpos)) deallocate(pekpos)
    
    ! if (measure_uncertainty) then
    !     allocate(eocpwr(num_xynod), eocden(num_xynod), asspwr(num_xynod,num_reload_comp+num_test_comp), &
    !         pekpwr(num_xynod,num_reload_comp+num_test_comp), pekass(num_xynod,num_reload_comp+num_test_comp), &
    !         pekden(num_xynod,num_reload_comp+num_test_comp), pekpos(2,num_xynod,num_reload_comp+num_test_comp), stat=ii)
    ! else
        allocate(eocpwr(num_xynod), eocden(num_xynod), asspwr(num_xynod,num_reload_comp), &
            pekpwr(num_xynod,num_reload_comp), pekfac(num_xynod,num_reload_comp), pekass(num_xynod,num_reload_comp), &
            pekden(num_xynod,num_reload_comp), pekpos(2,num_xynod,num_reload_comp), stat=ii)
    ! end if
    if(ii /= 0) then
        write(NOUT,100) 'failed to allocate arrays'
        call abort
    end if

    ! get EOC peak powers
    call get_peak_power(PWDINT, NOUT, eocpwr, peak_ass, peak_pwr, eocden, peknod, peak_pos, peak_den)
    eocpekpwr = peak_pwr
    eocpekdens = peak_den
    avgasspwr = 0._8
    DO ii = 1, num_potential_site
        IASS = fuel_pos(ii)
        avgasspwr = avgasspwr + eocpwr(IASS)
    END DO
    avgasspwr = avgasspwr / num_potential_site
    eocpekfac = eocpekpwr / avgasspwr
    ! edit EOC powers and find the second largest assembly/node powers 
    write(NOUT,'(/,20X,"*** Assembly powers at the end of cycle (EOC) ***")')
    if(num_hex_sector > 0) then
        write(NOUT,'("Index  IRING  IPOS   Power (W)    Peak density (W/cc)    Peak node",/)')
    else
        write(NOUT,'("Index    IX    IY    Power (W)    Peak density (W/cc)    Peak node",/)')
    end if
    peak_pwr2 = 0.0D0
    peak_den2 = 0.0D0
    IASS = 0
    do ii = 1, num_potential_site
        IASS = fuel_pos(ii)
        if(eocpwr(IASS) >= peak_pwr .and. IASS /= peak_ass) then
            peak_pwr2 = eocpwr(IASS)             ! second largest assembly power
            peak_ass2 = IASS
        end if
        if(eocden(IASS) >= peak_den .and. IASS /= peak_pos(1)) then 
            peak_den2 = eocden(IASS)             ! second largest power density in a different assembly
            peak_pos2(1) = IASS
            peak_pos2(2) = peknod(IASS)
        end if
        if(num_hex_sector > 0) then
            write(NOUT,'(I5,2I6,2ES15.5,7X,I6)') IASS, hexagon_rp(1:2,IASS), eocpwr(IASS), eocden(IASS), peknod(IASS)
        else
            call reload_ass2xy(IASS, num_xnode, IX, IY)
            write(NOUT,'(I5,2I6,2ES15.5,7X,I6)') IASS, IX, IY, eocpwr(IASS), eocden(IASS), peknod(IASS)
        end if
    end do
    write(NOUT,'(/,"Peak assembly power at EOC occurs in assembly ",I6, ", power (watts) =",ES13.5)') peak_ass, peak_pwr
    write(NOUT,'("Peak node-averaged power density at EOC occurs in assembly ",I6,", axial node ",I6,", density (w/cc) =",ES13.5,/)') &
    &   peak_pos(1:2), peak_den
    call reload_ass2xy(peak_ass, num_xnode, IX, IY)
    if(assembly_map(IX,IY) /= reload_type) write(NOUT,100) 'Warning: peak assembly power occurs at nonfuel assembly position'
    call reload_ass2xy(peak_pos(1), num_xnode, IX, IY)
    if(assembly_map(IX,IY) /= reload_type) write(NOUT,100) 'Warning: peak power density occurs at nonfuel assembly position'

    ! compute peak powers with reloaded assembly
    ! remember that now:
    !   peak_pwr / peak_pwr2  -  peak / second maximum assembly power at EOC
    !   peak_ass / peak_ass2  -  assembly position showing peak / second maximum assembly power
    !   peak_den / peak_den2  -  peak / second maximum (in a different assembly) node-averaged power density at EOC
    !   peak_pos / peak_pos2  -  (assembly, node) position showing peak power density / second maximum (in a different assembly)
    normfac = 1.0D0
    if(normalize_power) then
        write(NOUT,'("Renormalize assembly powers after refueling? Yes",/)')
    else
        write(NOUT,'("Renormalize assembly powers after refueling? No",/)')
    end if
    do iload = 1, num_reload_comp
        if(flux_data == 'RTFLUX') then
            do IY = 1, NINTJ 
                do IX = 1, NINTI
                    IASS = (IY-1)*NINTI + IX
                    asspwr(IASS,iload) = 0.0D0
                    pekpwr(IASS,iload) = 0.0D0
                    pekfac(IASS,iload) = 0.0D0
                    pekden(IASS,iload) = 0.0D0
                    local_peak_den = 0.0D0
                    if(assembly_map(IX, IY) /= reload_type) cycle
                    ! compute power for reloaded assembly under same flux
                    asspwr(IASS,iload) = eocpwr(IASS)
                    do IZ = bot_active_node, top_active_node
                        FDUM = sum(kerma%nheat(:,iload) * RTFLUX%FREG(IX,IY,IZ,:))
                        asspwr(IASS,iload) = asspwr(IASS,iload) + (FDUM - PWDINT%PWR(IX,IY,IZ)) * node_vol(IZ)
                        if(FDUM > local_peak_den) then
                            local_peak_den = FDUM
                            local_peak_nod = IZ
                        end if
                    end do
                    if(normalize_power) then
                        normfac = PWDINT%POWER / (PWDINT%POWER - eocpwr(IASS) + asspwr(IASS,iload))   ! renormalization factor
                        !write(NOUT,'("reloaded at",I6,"  renomalization factor =",ES13.5)') IASS, normfac
                    end if
                    asspwr(IASS,iload) = asspwr(IASS,iload) * normfac
                    local_peak_den = local_peak_den * normfac
                    if(IASS == peak_ass) then
                        new_peak_pwr = peak_pwr2 * normfac
                        new_peak_ass = peak_ass2
                    else
                        new_peak_pwr = peak_pwr * normfac
                        new_peak_ass = peak_ass
                    end if
                    if(IASS == peak_pos(1)) then
                        new_peak_den = peak_den2 * normfac 
                        new_peak_pos = peak_pos2
                    else
                        new_peak_den = peak_den * normfac
                        new_peak_pos = peak_pos
                    end if
                    if(asspwr(IASS,iload) > new_peak_pwr) then
                        new_peak_pwr = asspwr(IASS, iload)
                        pekass(IASS,iload) = IASS
                    else
                        pekass(IASS,iload) = new_peak_ass
                    end if
                    if(local_peak_den > new_peak_den) then
                        new_peak_den = local_peak_den
                        pekpos(:,IASS,iload) = (/IASS, local_peak_nod/)
                    else
                        pekpos(:,IASS,iload) = new_peak_pos
                    end if
                    pekpwr(IASS,iload) = new_peak_pwr       ! peak assembly-integrated power
                    pekfac(IASS,iload) = new_peak_pwr / avgasspwr
                    pekden(IASS,iload) = new_peak_den       ! peak node-averaged power density
                end do
            end do
        else if(flux_data == 'NHFLUX') then
            write(NOUT,100) 'Error: usage of NHFLUX not implemented yet'
            call abort
        end if
    end do

End Subroutine

! ----------------------------------------------------------------------
! get maximum increase and mean-square-root change in assembly power distribution from the
! reference distribution when a fresh assembly of type iload reloaded at position reload_pos
! Note: the power change is relative to the core-averaged assembly power
!       the maximum change consider only positive change (power increase)
Subroutine reload_power_change(iload, reload_pos, maxup, maxpos, rmsdif)
    use m_reload
    implicit none
    integer :: iload, reload_pos, maxpos
    real(8) :: maxup, rmsdif  ! max increase, root mean square difference
    ! local
    real(8) :: pwrdif, scale, inverse_avgpwr
    integer :: ii, IASS

    scale = (PWDINT%POWER - asspwr(reload_pos,iload)) / (PWDINT%POWER - eocpwr(reload_pos))
    inverse_avgpwr = num_potential_site*1.0d0 / PWDINT%POWER  
    ! it is fine to assume all power generated in fuel assemblies since this is used as a constant scale to quantify power difference
    maxup = 0.0d0
    rmsdif= 0.0d0
    do ii = 1, num_potential_site
        IASS = fuel_pos(ii)
        if(IASS == reload_pos) then
            pwrdif = inverse_avgpwr * (asspwr(reload_pos,iload) - refpwr(IASS))
        else
            pwrdif = inverse_avgpwr * (eocpwr(IASS)*scale - refpwr(IASS))
        end if
        if(pwrdif > maxup) then 
            maxup = pwrdif
            maxpos= IASS
        end if
        rmsdif = rmsdif + pwrdif*pwrdif
    end do
    rmsdif = sqrt(rmsdif/num_potential_site)
    
End Subroutine