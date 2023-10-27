Subroutine reload_distance(NOUT, num_candidate, candidate_pos)
    use m_reload
    implicit none
    integer :: NOUT, num_candidate
    integer :: candidate_pos(num_candidate)
    ! local
    integer :: iprev, ipos, IASS, JASS, IX, IY, JX, JY
    real(8) :: fdum, wgt(num_previous_site)
    real(8), external :: assembly_dist

    ! since cycle length and power is different, use burnup to determine weighting distance from each previous refueled assembly
    ! weight more the distance from assembly of smaller burnup
    fdum = 0.0d0
    do iprev = 1, num_previous_site
        IASS = previous_pos(iprev)
        wgt(iprev) = 1.0d0 / assbrn(IASS)
        fdum = fdum + wgt(iprev)
    end do
    wgt(:) = 1.0d0 / fdum * wgt(:)

    ! calculate weighted distance from previous refueled assemblies
    do ipos = 1, num_candidate
        IASS = candidate_pos(ipos)
        call reload_ass2xy(IASS, num_xnode, IX, IY)
        disperse_dist(IASS) = 0.0d0
        do iprev = 1, num_previous_site
            JASS = previous_pos(iprev)
            call reload_ass2xy(JASS, num_xnode, JX, JY)
            fdum = assembly_dist(IX,IY,JX,JY,num_hex_sector)
            disperse_dist(IASS) = disperse_dist(IASS) + wgt(iprev) * fdum
        end do
    end do

    !write(NOUT,'(20X,"*** Distance from recently refueled positions ***",/,&
    !&  "candidate position      distance (*assembly pitch)",/)')
    !write(NOUT,'((I4,ES13.5))') ((candidate_pos(ipos), disperse_dist(candidate_pos(ipos))), ipos=1,num_candidate)

End Subroutine