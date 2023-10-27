SUBROUTINE GetIJPOS(IRING,IPOS,IMESH,JMESH,NSECT,NCINTI,NCINTJ,NOUT)
    implicit none
    integer :: IRING,IPOS,IMESH,JMESH,NSECT,NCINTI,NCINTJ,NOUT
    !  returns (I,J) indices for a particular hexagon.
    !  IRING     intent(in)       RING NO. OF HEX.
    !  IPOS      intent(in)       POSITION NO. OF HEX IN RING.
    !  IMESH     intent(out)      1ST DIMENSION INDEX OF HEX.
    !  JMESH     intent(out)      2ND DIMENSION INDEX OF HEX.
    !  NSECT     intent(in)       NUMBER OF 60 DEGREE SECTORS IN PROBLEM DOMAIN. (= 1, 2 or 6)
    !  NCINTI, NCINTJ intent(in)     number of first and second dimension coarse mesh (node) intervals
    ! local
    integer :: ICENT,JCENT  ! center indices
    integer :: NDELI(6,2), NDELJ(6,2), NP, NS, IRM1, NPOS
    integer :: NTRIAG       ! rhombus type (hexagonal geometry, NTRIAG = 0 (1/3 or full core) or 1 (1/6 core)

    DATA NDELI / 0, -1, -1,  0,  1,  1,     -1, -1,  0,  1,  1,  0 /
    DATA NDELJ / 1,  0, -1, -1,  0,  1,      1,  0, -1, -1,  0,  1 / 
    !----------------------------------------------------------------------------------------------------------------------
    ! Two different rhombi:
    !                                                                         
    !     J (second dimension)                                                       J (second dimension)                                      
    !      \                                                                        /                                                   
    !       \                                                                      /                                                       
    !        ****************                                                     ****************                                      
    !         *              *                                                   *              *                                      
    !          *              *                                                 *              *                                     
    !           *              *                                               *              *                                    
    !            *              *                                             *              *                                   
    !             *              *                                           *              *                                  
    !              **************** ----> I (first dimension)               **************** ----> I (first dimension)       
    !                                                                                 
    !           corresponding to NDELI(:,1) and NDELJ(:,1)                  corresponding to NDELI(:,2) and NDELJ(:,2)                
    !------------------------------------------------------------------------------------------------------------------------ 
    100 format('[GetIJPOS]...',A,1X,I6)

    IF (NSECT == 2 .or. NSECT == 6) THEN
        NTRIAG = 1
    ELSE IF (NSECT == 1) THEN 
        NTRIAG = 2
    ELSE
        write(NOUT,100) 'Error: wrong NSECT input to GetIJPOS, in hex node, NSECT = 1, 2 or 6'
    END IF

    ! DETERMINE (ICENT,JCENT), THE INDICES OF THE CENTRAL HEXAGON.
    ICENT = (NCINTI+1)/2
    JCENT = ICENT
    IF (NSECT .NE. 6) THEN
        ICENT = 1
        JCENT = 1
    END IF

    !     USING INDICES OF CENTRAL HEX, SET INDICES OF FIRST HEX IN RING IRING.
    IMESH = ICENT + IRING - 1
    JMESH = JCENT

    IF (IRING .GT. 1) THEN
        ! MOVE AROUND RING TO FIND HEX IN POSITION IPOS.
        IRM1 = IRING - 1
        NPOS = 1
        !DO NS=1,NSECT
        DO NS=1,6
            DO NP=1,IRM1
                IF (NPOS .NE. IPOS) THEN
                    NPOS = NPOS + 1
                    ! Because NTRIAG = 1 / 2 only (hexagonal nodal case, when the nodal mapping is needed)
                    ! Only two dimensional NDELI and NDELJ are needed.
                    ! They are used to express the (I,J) indices' change in 6 sections along each ring. 
                    IMESH = IMESH + NDELI(NS,NTRIAG)
                    JMESH = JMESH + NDELJ(NS,NTRIAG)
                END IF
            END DO
        END DO
    END IF
    !!!!!   IF HEX LIES OUTSIDE GEODST DOMAIN, SET IMESH (JMESH) = -1.
    IF (IMESH .LT. 1 .OR. IMESH .GT. NCINTI) IMESH = -1
    IF (JMESH .LT. 1 .OR. JMESH .GT. NCINTJ) JMESH = -1

END SUBROUTINE

! ----------------------------------------------------------------------
SUBROUTINE GetHexMap(ITRMAP,NHEX,NSECT,NCINTI,NCINTJ,NOUT)
    !  CONSTRUCTS THE TRANSFORMATION MAP BETWEEN THE GEODST AND NODAL MESH CELL ORDERINGS.
    !  ITRMAP   intent(out)         ITRMAP(IHEX) = IJ, WHERE IJ AND IHEX ARE THE
    !                               NUMBERS OF A GIVEN HEXAGON IN GEODST AND NODAL
    !                               ORDERING, RESPECTIVELY.  IJ IS SET TO  0
    !                               FOR HEXAGONS LYING OUTSIDE THE GEODST DOMAIN.
    !  NHEX     intent(out)         Total number of nodes in problem domain
    implicit none
    integer :: NSECT,NCINTI,NCINTJ,NOUT
    integer :: ITRMAP(NCINTI*NCINTJ)
    ! local
    integer :: NHEX, IPOS, IRING, I, J, IJ, NH, NRINGT

    NRINGT = NCINTI
    IF (NSECT == 6) NRINGT = (NCINTI+1)/2

    !     CONSTRUCT ITRMAP IN HEXAGONAL GEOMETRY.
    NHEX = 0
    DO IRING=1,NRINGT
        NH = NSECT*(IRING-1)
        NH = MAX(NH,1)
        DO IPOS=1,NH
            NHEX = NHEX + 1
            CALL GetIJPOS(IRING,IPOS,I,J,NSECT,NCINTI,NCINTJ,NOUT)
            IJ = (J-1)*NCINTI + I
            IF (I < 0 .OR. J < 0) IJ = 0
            ITRMAP(NHEX) = IJ
        END DO
    END DO

END SUBROUTINE

! ----------------------------------------------------------------------
Subroutine MapIJ2RP(ITRMAP2, NSECT, NCINTI, NCINTJ, NOUT)
    ! get Ring and Position numbers for given GEODST index
    integer :: NSECT, NCINTI, NCINTJ, NOUT
    integer :: ITRMAP2(2,NCINTI*NCINTJ)     
    ! ITRMAP2(1, IJ) = ring number,  ITRMAP2(2, IJ) = position number
    ! local
    integer :: IR, IP, NR, NH, I, J, IJ

    NR = NCINTI 
    if(NSECT == 6) NR = (NR+1)/2

    ITRMAP2 = 0
    do IR = 1, NR 
        NH = NSECT * (IR-1)
        NH = max(NH,1)
        do IP = 1, NH
            call GetIJPOS(IR, IP, I, J, NSECT, NCINTI, NCINTJ, NOUT)
            if (I < 0 .or. J < 0) cycle
            IJ = (J-1)*NCINTI + I
            ITRMAP2(1,IJ) = IR
            ITRMAP2(2,IJ) = IP
        end do
    end do

End Subroutine

! ----------------------------------------------------------------------
! compute distance between two assemblies, in unit of assembly pitch
! considering hexagonal core or Cartesian core with square lattice
Function assembly_dist(IX,IY,JX,JY,NSECT) result(dist)
    implicit none 
    integer :: IX, IY, JX, JY, NSECT
    real(8) :: dx, dy, dist, sqrt3
    data sqrt3 /1.7320508d0/

    dx = IX - JX
    dy = IY - JY
    if(NSECT == 0) then   ! Cartesian (square)
        dist = sqrt(dx*dx + dy*dy)
    else if(NSECT == 1) then 
        dx = dx + dy*0.5d0
        dy = dy * sqrt3 * 0.5d0
        dist = sqrt(dx*dx + dy*dy)
    else if(NSECT == 2 .or. NSECT == 6) then
        dx = dx - dy*0.5d0
        dy = dy * sqrt3 * 0.5d0
        dist = sqrt(dx*dx + dy*dy)
    else
        write(0,'("[assembly_dist]...Error: invalid parameter NSECT",I4)') NSECT
        call abort
    end if
End Function
