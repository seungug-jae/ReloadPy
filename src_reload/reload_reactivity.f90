! estimate reactivity change due to fuel assembly reloading
! it is time consuming even with FOP, consider to improve it by not perturbing all fuel assemblies but only candidates
Subroutine reload_reactivity(NOUT, num_candidate, candidate_pos)
    use m_reload
    use m_utility
    use NE_Kind, only: NE_Kind_GetFreeLogicalUnit
    implicit none
    integer :: NOUT, num_candidate, candidate_pos(num_candidate)
    ! local
    integer :: ii, jj, kk, IASS
    integer :: NIN, MCARD, ncards(99), freeform, HMN
    integer :: chdir    ! intrinsic function
    character(len=8) :: ANAME
    logical :: forward_avaiable = .false.
    logical :: adjoint_avaiable = .false.
    logical :: adjoint_initial = .false.
    real(8) :: fdum, fdum2

    100 format('[reload_reactivity]...',A,1X,I6)

    if(use_persent) then
        ! check VARIANT option
        inquire(file='user.NHFLUX', exist=forward_avaiable)
        inquire(file='user.NAFLUX', exist=adjoint_avaiable)
        if(.not. forward_avaiable) then
            NIN = NE_Kind_GetFreeLogicalUnit()
            open(NIN,file='ADIF3D',status='old',form='formatted')
            call scan_arc_file_header(NIN,MCARD,ncards,freeform,ANAME)
            if(freeform /= 1) then
                write(NOUT,100) 'Error: when PERSENT is invoked, A.DIF3D input has to be in free format'
                call abort 
            end if
            if(MCARD >= 12 .and. ncards(12) > 0) then 
                ii = sum(ncards(1:11))
                call SkipLines(NIN, ii)
                read(NIN,*) ii, ii, HMN
                HMN = abs(HMN)
                if(HMN > 9999) HMN = modulo(HMN,10000)
                if(HMN > 99) then
                    ii = HMN / 100              ! M = N is guaranteed 
                else
                    ii = HMN / 10
                end if
                if(ii == 1) then
                    forward_avaiable = .true.     ! see below, if NHFLUX at EOC is saved, we can fix this 
                    call system_command(NOUT,'ln -s  NHFLX0  user.NHFLUX')
                end if
            end if
            close(NIN)
        end if

        ! set up persent run environment
        if(adjoint_avaiable) then
            call system_command(NOUT, 'rm '//persent_dir//'/*')    ! not first run
        else
            call system_command(NOUT, 'mkdir '//persent_dir)
            inquire(file='init.NAFLUX', exist=adjoint_initial)
            if(adjoint_initial) call system_command(NOUT, 'cp init.NAFLUX '//persent_dir//'/NAFLUX')   ! used as initial guess for adjoint calculation
        end if
        call system_command(NOUT, 'cp ISOTXS '//persent_dir//'/ISOTXS')         ! for DIF3D null run
        call system_command(NOUT, 'cp ISOTXS '//persent_dir//'/user.ISOTXS')    ! for PERSENT
        call system_command(NOUT, 'cp GEODST '//persent_dir//'/GEODST')
        call system_command(NOUT, 'cp NDXSRF '//persent_dir//'/NDXSRF')
        call system_command(NOUT, 'cp LABELS '//persent_dir//'/LABELS')
        call system_command(NOUT, 'cp ZNATDN '//persent_dir//'/ZNATDN')

        ! try to utilize REBUS output, but it seems REBUS has problem with saving NHFLUX at EOC. Only NHFLX0 is updated at EOC
        if(forward_avaiable) then
            if(adjoint_avaiable) then
                call reload_persent_input('persent.inp',2, num_candidate, candidate_pos)
            else
                call reload_persent_input('persent.inp',1, num_candidate, candidate_pos)
            end if
        else
            call reload_persent_input('persent.inp',0, num_candidate, candidate_pos)   ! let PERSENT generate NHFLUX and NAFLUX
        end if
        ! make dif3d input using interface files
        if(adjoint_initial) then
            call reload_dif3d_input(NOUT,'reload.dif3d',2)
        else
            call reload_dif3d_input(NOUT,'reload.dif3d',0)
        end if
        call system_command(NOUT, 'cp reload.dif3d '//persent_dir//'/dif3d.inp')
        call system_command(NOUT, 'cp persent.inp '//persent_dir//'/persent.inp')
        
        ! move in temporary working directory
        ii = chdir(persent_dir)
        if(ii /= 0) then
            write(NOUT,100) 'Error: failed to change working directory to '//persent_dir
            call abort 
        end if

        ! link existing flux data (avoid copy of large files)
        if(forward_avaiable) call system_command(NOUT, 'ln -s ../user.NHFLUX .')
        if(adjoint_avaiable) call system_command(NOUT, 'ln -s ../user.NAFLUX .')
        
        ! run persent
        call CPU_TIME(fdum)
        call system_command(NOUT, '../persent.x < persent.inp > persent.out')
        call CPU_TIME(fdum2)
        persent_time(2) = persent_time(2) + fdum2 - fdum
        ! grab reactivity perturbation results for each fuel assembly from persent output
        call CPU_TIME(fdum)
        call reload_persent_output(NOUT, 'persent.out', eoc_keff)
        call CPU_TIME(fdum2)
        persent_time(3) = persent_time(3) + fdum2 - fdum

        ! save NHFLUX and NAFLUX in case additional PERSENT run is needed
        if(.not. forward_avaiable) call system_command(NOUT, 'mv base.NHFLUX ../user.NHFLUX')
        if(.not. adjoint_avaiable) call system_command(NOUT, 'mv base.NAFLUX ../user.NAFLUX')
        
        ! compute assembly worth 
        if(consider_worth) then
            if(.not. allocated(asswth)) then
                allocate(asswth(num_xynod, num_reload_comp))
                asswth(:,:) = 0.0d0
            end if
            jj = num_reload_comp + 1
            do ii = 1, num_reload_comp
                do kk = 1, num_candidate
                    IASS = candidate_pos(kk)
                    asswth(IASS,ii) = deltak(IASS,ii) - deltak(IASS,jj)
                end do
            end do
        end if

        ! move back to REBUS working directory
        ii = chdir('..')
        if(ii /= 0) then
            write(NOUT,100) 'Error: failed to move back to REBUS working directory'
            call abort 
        end if
    end if    

End Subroutine
