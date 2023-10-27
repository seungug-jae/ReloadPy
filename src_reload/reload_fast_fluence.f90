! compute accumulated nodal average fast fluence in one cycle time interval as
! fluence = 1/2(boc_fast_flux + eoc_fast_flux) * cycle_length
! after each cycle, node averaged fast fluence will be saved into binary file reload.fluence
! format: 
!   1D    NINTI, NINTJ, NINTK, irradiation_time (day)
!   2D    fluence(NINTI,NINTJ,NINTK)
Subroutine reload_fast_fluence(NOUT, NINTI, NINTJ, NINTK, NGROUP, bocflx, eocflx)
    use m_reload    !, only: fluence, assembly_map, reload_type, cycle_length, grpstr, min_fast_energy, num_hex_sector, hexagon_rp, asscyc
    implicit none
    integer :: NOUT, NINTI, NINTJ, NINTK, NGROUP
    real(8) :: bocflx(NINTI,NINTJ,NINTK,NGROUP), eocflx(NINTI,NINTJ,NINTK,NGROUP)
    ! local
    integer, parameter :: day = 86400    ! 3600*24 = 86400 sec
    real(8) :: old_fluence(NINTI,NINTJ,NINTK), time
    integer :: gend, gg, NIN, IOS, NI, NJ, NK, IASS, peak_nod, peak_ass
    real(8) :: weight, peak_fluence, global_peak
    logical :: yes
    data NIN /100/

    100 format('[reload_fast_fluence]...',A,1X,I6)

    if(min_fast_energy >= grpstr(1)) then
        write(NOUT,100) 'Error: cut-off energy of fast neutron is greater than group 1 upper boundary'
        call abort
    end if
    do gg = 2, NGROUP+1
        if(grpstr(gg) <= min_fast_energy) then
            gend = gg - 1
            exit
        end if 
    end do
    if(gg > NGROUP+1) gend = NGROUP
    weight = log(grpstr(gend)/min_fast_energy) / log(grpstr(gend)/grpstr(gend+1))

    if(allocated(fluence)) deallocate(fluence)
    allocate(fluence(NINTI,NINTJ,NINTK))
    fluence(:,:,:) = 0.0D0
    do gg = 1, gend-1
        fluence = fluence + bocflx(:,:,:,gg) + eocflx(:,:,:,gg)
    end do
    fluence = fluence + weight * (bocflx(:,:,:,gend) + eocflx(:,:,:,gend))
    fluence = 0.5D0 * fluence * cycle_length * day

    inquire(file='reload.fluence', exist=yes)
    if(yes) then
        open(NIN,file='reload.fluence',status='old',form='unformatted',iostat=IOS)
        if(IOS /= 0) then
            write(NOUT,100) 'Error: cannot open existing file reload.fluence'
            call abort
        end if
        read(NIN,iostat=IOS) NI, NJ, NK, time
        if(IOS /= 0) then
            write(NOUT,100) 'Error: I/O error when reading file reload.fluence record', 1
            call abort
        end if
        if(NI /= NINTI .or. NJ /= NINTJ .or. NK /= NINTK) then
            write(NOUT,100) 'Error: the existing fluence data has inconsistent dimensions'
            call abort 
        end if
        read(NIN,iostat=IOS) old_fluence
        if(IOS /= 0) then
            write(NOUT,100) 'Error: I/O error when reading file reload.fluence record', 2
            call abort
        end if
        close(NIN)
        fluence = fluence + old_fluence
    end if

    ! save updated fast fluence
    open(NIN,file='reload.fluence',status='replace',form='unformatted',iostat=IOS)
    if(IOS /= 0) then
        write(NOUT,100) 'Error: cannot create file reload.fluence'
        call abort
    end if
    write(NIN) NINTI, NINTJ, NINTK, time+cycle_length
    write(NIN) fluence
    close(NIN)

    ! edit peak fast fluence
    write(NOUT,'(/,20X,"*** Peak node averaged fast fluence (E > ",F6.2," keV) in fuel assembly ***")') min_fast_energy*1e-3
    if(num_hex_sector > 0) then
        write(NOUT,'("Index  IRING  IPOS   Peak fluence (n/cm^2)   Burned cycles   Peak node",/)')
    else
        write(NOUT,'("Index    IX    IY    Peak fluence (n/cm^2)   Burned cycles   Peak node",/)')
    end if
    IASS = 0
    global_peak = 0.0D0
    do NJ = 1, NINTJ
        do NI = 1, NINTI
            IASS = IASS + 1
            if(assembly_map(NI,NJ) /= reload_type) cycle
            peak_fluence = 0.0D0
            do NK = 1, NINTK
                if(fluence(NI,NJ,NK) > peak_fluence) then
                    peak_fluence = fluence(NI,NJ,NK)
                    peak_nod = NK
                end if
            end do
            if(num_hex_sector > 0) then
                write(NOUT,'(I5,2I6,6X,ES15.5,2(9X,I6))') IASS, hexagon_rp(1:2,IASS), peak_fluence, asscyc(IASS), peak_nod
            else
                write(NOUT,'(I5,2I6,6X,ES15.5,2(9X,I6))') IASS, NI, NJ, peak_fluence, asscyc(IASS), peak_nod
            end if
            if(peak_fluence > global_peak) then
                global_peak = peak_fluence
                peak_ass = IASS
            end if
        end do
    end do
    write(NOUT,'(/,"Peak node averaged fast fluence occurs in assembly ",I6,", fluence (n/cm^2) =",ES13.5)') peak_ass, global_peak

End subroutine
