Subroutine get_peak_power(PWDINT, NOUT, pwr, peak_ass, peak_pwr, pwrden, peak_nod, peak_pos, peak_den)
    ! calculate assembly powers, peak assembly power and position, determine peak node power density in each assembly
    ! and position, and the global peak node-averaged power density using PWDINT of DIF3D/REBUS output
    use m_reload, only: PWDINT_DATA, num_xnode, num_ynode, num_xynod, node_vol, assembly_map
    implicit none
    type(PWDINT_DATA) :: PWDINT
    integer :: NOUT, peak_ass, peak_pos(2), peak_nod(num_xynod)
    real(8) :: pwr(num_xynod), peak_pwr     ! assembly power, peak assembly power
    real(8) :: pwrden(num_xynod), peak_den  ! peak power density (node level for now)
    ! local
    integer :: IX, IY, IZ, IASS
    
    100 format('[get_peak_power]...',A,1X,I6)

    if(.not. PWDINT%defined) then
        write(NOUT,100) 'Error: PWDINT data not imported'
        call abort
    end if
    if(.not. allocated(node_vol)) then
        write(NOUT,100) 'Error: node volumes (node_vol) not computed'
        call abort
    end if

    peak_pwr = 0.0D0
    peak_den = 0.0D0
    IASS = 0
    do IY = 1, num_ynode
        do IX = 1, num_xnode
            IASS = IASS + 1
            if(assembly_map(IX,IY) == 0) then
                pwr(IASS) = 0.0D0 
                pwrden(IASS) = 0.0D0
                peak_nod(IASS) = 0
                cycle
            end if
            pwr(IASS) = sum(PWDINT%PWR(IX,IY,:) * node_vol(:))
            if(pwr(IASS) >= peak_pwr) then
                peak_pwr = pwr(IASS)
                peak_ass = IASS
            end if
            pwrden(IASS) = 0.0D0
            do IZ = 1, PWDINT%NINTK
                if(PWDINT%PWR(IX,IY,IZ) > pwrden(IASS)) then
                    pwrden(IASS) = PWDINT%PWR(IX,IY,IZ)
                    peak_nod(IASS) = IZ
                end if
            end do
            if(pwrden(IASS) > peak_den) then
                peak_den = pwrden(IASS)
                peak_pos(1) = IASS
                peak_pos(2) = peak_nod(IASS)
            end if
        end do
    end do

End Subroutine