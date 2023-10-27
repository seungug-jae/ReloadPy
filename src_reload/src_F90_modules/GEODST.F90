!---------------------------------------------------------------------------------------------------------------------------------
!  Using the historical format for the GEODST file, import/export/print the data from the binary file
!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_ASSIGNPRINTINFO           !  Sets the output unit for the module
!  GEODST_DEFINE                    !  Defines a GEODST_DATA type based upon user input
!  GEODST_Copy                      !  Makes a copy of an existing GEODST data structure
!  GEODST_NRASS1                    !  Converts the coarse mesh mapping (NRASS=0) to fine mesh mapping (NRASS1) in GEODST_DATA
!  GEODST_AddMeshes                 !  This routine adds the user specified meshes to an existing domain
!  GEODST_LocatePoint               !  This routine identifies the GEODST I,J,K position and local coordinates of a given point within that mesh
!  GEODST_ComputeCentroids          !  This routine computes and returns the centroids of each region in the GEODST mesh
!  GEODST_ComputeAverageThickness   !  This routine computes and returns the average thickness of each region in the GEODST mesh with respect to the coordinate direction
!  GEODST_MeshVolumes               !  This routine computes and returns the mesh by mesh volume in the GEODST mesh
!  GEODST_RegionsOnBoundary         !  This routine computes and returns a vector that defines which regions touch which boundaries and roughly at which position
!  GEODST_NODAL_factors             !  This routine computes and returns the mesh by mesh NODAL flux factors
!  GEODST_VOID                      !  Provides a path to void a GEODST data structure
!  GEODST_PRINT                     !  Prints the entire GEODST data structure
!-----------------------------------
!  GEODST_ASSIGNPRINTINFO(MODULE_OUTPUT_UNIT)
!  GEODST_DEFINE(User_GEODST,USER_IGOM,USER_NZONE,USER_NREG,USER_NCINTI,USER_NCINTJ,USER_NCINTK,        &
!                     USER_NINTI,USER_NINTJ,USER_NINTK,USER_NBS,USER_NBCS,USER_NIBCS,USER_NZWBB,USER_NRASS   )
!  GEODST_Copy(Source_GEODST,Destination_GEODST)
!  GEODST_NRASS1(User_GEODST)
!  GEODST_AddMeshes(User_GEODST,IJK,NumBoundaries,MeshBoundaries)
!  GEODST_LocatePoint(User_GEODST,X,Y,Z,oINode,oJNode,oKNode,o_Reg,oX,oY,oZ)
!  GEODST_ComputeCentroids(User_GEODST,Centroids)
!  GEODST_ComputeAverageThickness(User_GEODST,Thickness)
!  GEODST_MeshVolumes(User_GEODST,MeshVolumes)
!  GEODST_RegionsOnBoundary(User_GEODST,Coordinates)
!  GEODST_NODAL_factors(User_GEODST,MeshFluxFactors)
!  GEODST_VOID(User_GEODST)
!  GEODST_PRINT(User_GEODST)
!---------------------------------------------------------------------------------------------------------------------------------

MODULE GEODST_io 
#include "DIF3D_Types.h"
IMPLICIT NONE

TYPE GEODST_DATA
   DIF3D_Log :: DEFINED=.FALSE.                      ! Defines whether the data structure is defined
   DIF3D_FileNames       :: FILENAME = 'UNDEFINED'    ! The originating file name
   ! CARD TYPE 0   FILE IDENTIFICATION 
   DIF3D_Char        :: HNAME = 'UNDEFINE'
   DIF3D_Char        :: HUSE(2) = 'UNDEFINE'
   GEODST_Int :: IVERS = 0
   ! CARD TYPE 1           ---FILE SPECIFICATIONS---
   GEODST_Int :: IGOM = -1  ! GEOMETRY TYPE 
                                    !       0- POINT (FUNDAMENTAL MODE),  1- SLAB  ,  2- CYLINDER  ,  3- SPHERE               ,
                                    !       6- X-Y                     ,  7- R-Z   ,  8- THETA-R   ,  9- UNIFORM TRIANGULAR   ,
                                    !      10- HEXAGONAL (1 MESH POINT IN EACH HEXAGONAL ELEMENT)  , 11- R-THETA              ,
                                    !      12- R-THETA-Z,              , 13- R-THETA-ALPHA         , 14- X-Y-Z                ,
                                    !      15- THETA-R-Z,              , 16- THETA-R-ALPHA,        , 17- UNIFORM TRIANGULAR-Z ,
                                    !      18- HEXAGON-Z   (MESH POINTS AS IN 10 ABOVE) 
   GEODST_Int :: NZONE      ! NUMBER OF ZONES (EACH HOMOGENEOUS IN NEUTRONICS PROBLEM - A ZONE CONTAINS ONE OR MORE REGIONS)
   GEODST_Int :: NREG       ! NUMBER OF REGIONS
   GEODST_Int :: NZCL       ! NUMBER OF ZONE CLASSIFICATIONS (EDIT PURPOSES)
   GEODST_Int :: NCINTI     ! NUMBER OF FIRST DIMENSION COARSE MESH INTERVALS
   GEODST_Int :: NCINTJ     ! NUMBER OF SECOND DIMENSION COARSE MESH INTERVALS. NCINTJ.EQ.1 FOR ONE DIMENSIONAL CASE.
   GEODST_Int :: NCINTK     ! NUMBER OF THIRD DIMENSION COARSE MESH INTERVALS NCINTK.EQ.1 FOR ONE AND TWO DIMENSIONAL CASES.
   GEODST_Int :: NINTI      ! NUMBER OF FIRST DIMENSION FINE MESH INTERVALS
   GEODST_Int :: NINTJ      ! NUMBER OF SECOND DIMENSION FINE MESH INTERVALS NINTJ.EQ.1 FOR ONE DIMENSIONAL CASE.
   GEODST_Int :: NINTK      ! NUMBER OF THIRD DIMENSION FINE MESH INTERVALS NINTK.EQ.1 FOR ONE AND TWO DIMENSION CASES.
   GEODST_Int :: IMB1       ! FIRST BOUNDARY ON FIRST DIMENSION
   GEODST_Int :: IMB2       ! LAST BOUNDARY ON FIRST DIMENSION
   GEODST_Int :: JMB1       ! FIRST BOUNDARY ON SECOND DIMENSION
   GEODST_Int :: JMB2       ! LAST BOUNDARY ON SECOND DIMENSION
   GEODST_Int :: KMB1       ! FIRST BOUNDARY ON THIRD DIMENSION
   GEODST_Int :: KMB2       ! LAST BOUNDARY ON THIRD DIMENSION
                                    !    BOUNDARY DEFINITIONS 0 - ZERO FLUX (DIFFUSION)
                                    !                         1 - REFLECTED
                                    !                         2 - EXTRAPOLATED (DIFFUSION - DEL PHI/PHI = -C/D WHERE C IS GIVEN 
                                    !                             AS BNDC BELOW AND D IS THE GROUP DIFFUSION CONSTANT, TRANSPORT - NO RETURN).
                                    !                         3 - REPEATING (PERIODIC) WITH OPPOSITE FACE
                                    !                         4 - REPEATING (PERIODIC) WITH NEXT ADJACENT FACE.
                                    !                         5 - INVERTED REPEATING ALONG THIS FACE. (180 DEGREE ROTATION)
                                    !                         6 - ISOTROPIC RETURN (TRANSPORT)
                                    !         NOTE FOR REPEATING CONDITIONS (3,4,5) - LET I1 DENOTE FIRST BOUNDARY ON FIRST DIMENSION, I2 THE SECOND BOUNDARY ON THE FIRST DIMENSION, 
                                    !         J1 THE FIRST BOUNDARY ON THE SECOND DIMENSION, ETC. THEN THESE REPEATING BOUNDARY CONDITIONS ONLY APPLY TO BOUNDARIES I1,I2,J1, AND J2.
                                    !         GOING IN ORDER OF I1,J1,I2,J2, THE FIRST BOUNDARY WHICH IS INVOLVED CARRIES THE DESIGNATOR DEFINING THE REPEATING CONDITION.
   GEODST_Int :: NBS        ! NUMBER OF BUCKLING SPECIFICATIONS
                                    !      0 - NONE
                                    !      1 - SINGLE VALUE APPLIES EVERYWHERE.EQ.NZONE-  ZONE DEPENDENT M*NZONE - DATA IS GIVEN OVER ALL ZONES FOR
                                    !          THE FIRST ENERGY GROUP, THEN FOR THE NEXT GROUP, TO END OF LIST.  
                                    !          IF M.LT.NGROUP THEN THE M-TH GROUP DATA APPLIES TO ALL ADDITIONAL GROUPS. (2.LE.M.LE.NGROUP)
   GEODST_Int :: NBCS       ! NUMBER OF CONSTANTS FOR EXTERNAL BOUNDARIES
                                    !      0   - NONE
                                    !      1   - SINGLE VALUE USED EVERYWHERE
                                    !      6   - INDIVIDUAL VALUE GIVEN FOR EACH EXTERNAL BOUNDARY. 
                                    !            THE ORDERING OF THE VALUES IS THE SAME AS THE ORDERING OF THE BOUNDARY CONDITIONS.
                                    !      6*M - SIX VALUES GIVEN FOR FIRST ENERGY GROUP (ORDERED AS DESCRIBED ABOVE),
                                    !            THEN 6 FOR THE NEXT GROUP, TO END OF LIST.  (2.LE.M.LE.NGROUP).
                                    !            IF M.LT.NGROUP THEN THE M-TH GROUP DATA APPLIES TO ALL REMAINING GROUPS.
   GEODST_Int :: NIBCS      ! NUMBER OF CONSTANTS FOR INTERNAL BOUNDARIES
                                    !      0 - NONE
                                    !      1 - SINGLE VALUE USED EVERYWHERE 1.LT.N.LE.NGROUP 
                                    !          -VALUES ARE GIVEN BY ENERGY GROUP WITH NON-BLACK CONDITION INDICATED BY ZERO ENTRY- 
                                    !          LAST VALUE APPLIES TO ADDITIONAL GROUPS
                                    !          ***ANL MOD***     NZWBB*NGROUP - ZONE- AND GROUP-DEPENDENT DATA. THE N-TH SET OF NGROUP GROUP-DEPENDENT
                                    !          ***ANL MOD***     DATA CONSTANTS IS ASSIGNED TO THE BLACK ABSORBER ZONE NZHBB(N).
   GEODST_Int :: NZWBB      ! NUMBER OF ZONES WHICH ARE BLACK ABSORBERS
   GEODST_Int :: NTRIAG     ! TRIANGULAR/HEXAGONAL GEOMETRY OPTION
                                    !      0 - REGION OF SOLUTION IS A RHOMBUS IN WHICH THE 1ST AND 2ND DIMENSION AXES INTERSECT AT AN ANGLE OF 120 DEGREES.
                                    !      1 - REGION OF SOLUTION IS A RHOMBUS IN WHICH THE 1ST AND 2ND DIMENSION AXES INTERSECT AT AN ANGLE OF 60 DEGREES.
                                    !      2 - REGION OF SOLUTION IS A RECTANGLE.  THE BOUNDARIES I1 AND I2 BISECT MESH TRIANGLES. SEE NTHPT BELOW. (IGOM=9,17 ONLY)
                                    !      3 - REGION OF SOLUTION IS AN EQUILATERAL, 60 DEGREE TRIANGLE.  (IGOM=9,17 ONLY)
                                    !      4 - REGION OF SOLUTION IS A 30-60 DEGREE RIGHT TRIANGLE IN WHICH THE 1ST AND 2ND 
                                    !          DIMENSION AXES INTERSECT AT THE 30 DEGREE ANGLE.  (IGOM=9,17 ONLY)
                                    !      5 - REGION OF SOLUTION IS A RHOMBUS IN WHICH THE 1ST AND 2ND DIMENSION AXES INTERSECT 
                                    !          AT AN ANGLE OF 30 DEGREES. (IGOM=9,17 ONLY)
   GEODST_Int :: NRASS      ! REGION ASSIGNMENTS
                                    !      0- TO COARSE MESH
                                    !      1- TO FINE MESH
   GEODST_Int :: NTHPT      ! ORIENTATION OF FIRST FINE MESH INTERVAL IN TRIANGULAR GEOMETRIES.  NTRIAG=2 ONLY.
                                    !      1- TRIANGLE(1,1) POINTS AWAY FROM FIRST DIMENSION AXIS, I.E., NO INTERNAL MESH LINE INTERSECTS THE ORIGIN.
                                    !      2- TRIANGLE(1,1) POINTS TOWARD THE FIRST DIMENSION AXIS, I.E., AN INTERNAL MESH LINE INTERSECTS THE ORIGIN.
   GEODST_Int :: NGOP(4)    ! RESERVED

   ! CARD TYPES 2-4: UNITS ARE CM FOR LINEAR DIMENSIONS AND RADIANS FOR ANGULAR DIMENSIONS
   !  FOR UNIFORM-TRIANGULAR-MESH GEOMETRY (IGOM = 9) THE LENGTH (L) OF THE SIDE OF A MESH TRIANGLE MUST BE GIVEN BY THE EXPRESSION: L = 2.*(XMESH(2)-XMESH(1))/IFINTS(1) . 
   !  FOR UNIFORM-HEXAGONAL-MESH GEOMETRY (IGOM = 10) THE FLAT-TO-FLAT DISTANCE (FTF) ACROSS A MESH HEXAGON MUST BE GIVEN BY THE EXPRESSION FTF = (XMESH(2)-XMESH(1))/IFINTS(1)
   !  FOR UNIFORM-TRIANGULAR-MESH GEOMETRY (IGOM = 17) THE LENGTH (L) OF THE SIDE OF A MESH TRIANGLE MUST BE GIVEN BY THE EXPRESSION: L = 2.*(XMESH(2)-XMESH(1))/IFINTS(1) .
   !  FOR UNIFORM-HEXAGONAL-MESH GEOMETRY (IGOM = 18) THE FLAT-TO-FLAT DISTANCE (FTF) ACROSS A MESH HEXAGON MUST BE GIVEN BY THE EXPRESSION FTF = (XMESH(2)-XMESH(1))/IFINTS(1)

   ! CARD TYPE 2           ---ONE DIMENSIONAL COARSE MESH INTERVAL BOUNDARIES AND FINE MESH INTERVALS (2D RECORD)---
   ! PRESENT IF IGOM.GT.0 AND IGOM.LE.3
   GEODST_Real, DIMENSION(:), POINTER :: XMESH   ! NCINTI+1           ! COARSE MESH BOUNDARIES, FIRST DIMENSION
   GEODST_Int,  DIMENSION(:), POINTER :: IFINTS            ! NCINTI             ! NUMBER OF EQUALLY SPACED FINE MESH INTERVALS PER COARSE MESH INTERVAL, FIRST DIMENSION.

   ! CARD TYPE 3           ---TWO DIMENSIONAL COARSE MESH INTERVAL BOUNDARIES AND FINE MESH INTERVALS (3D RECORD)---
   ! PRESENT IF IGOM.GE.6 AND IGOM.LE.11
   GEODST_Real, DIMENSION(:), POINTER :: YMESH   ! NCINTJ+1           ! COARSE MESH BOUNDARIES, SECOND DIMENSION
   GEODST_Int,  DIMENSION(:), POINTER :: JFINTS            ! NCINTJ             ! NUMBER OF EQUALLY SPACED FINE MESH INTERVALS PER COARSE MESH INTERVAL, SECOND DIMENSION.

   ! CARD TYPE 4           ---THREE DIMENSIONAL COARSE MESH INTERVAL BOUNDARIES AND FINE MESH INTERVALS (4D RECORD)---
   ! PRESENT IF IGOM.GE.12
   GEODST_Real, DIMENSION(:), POINTER :: ZMESH   ! NCINTK+1           ! COARSE MESH BOUNDARIES, THIRD DIMENSION
   GEODST_Int,  DIMENSION(:), POINTER :: KFINTS            ! NCINTK             ! NUMBER OF EQUALLY SPACED FINE MESH INTERVALS PER COARSE MESH INTERVAL, THIRD DIMENSION.

   ! CARD TYPE 5           ---GEOMETRY DATA (5D RECORD)---
   ! PRESENT IF IGOM.GT.0 OR NBS.GT.0                           -
   GEODST_Real_5D, DIMENSION(:), POINTER :: VOLR   ! NREG               ! REGION VOLUMES (CC)
   GEODST_Real_5D, DIMENSION(:), POINTER :: BSQ    ! NBS                ! BUCKLING (B**2) VALUES (CM**-2)
   GEODST_Real_5D, DIMENSION(:), POINTER :: BNDC   ! NBCS               ! BOUNDARY CONSTANTS (DEL PHI/PHI =-C/D)
   GEODST_Real_5D, DIMENSION(:), POINTER :: BNCI   ! NIBCS              ! INTERNAL BLACK BOUNDARY CONSTANTS
   GEODST_Int , POINTER :: NZHBB(:)  ! NZWBB              ! ZONE NUMBERS WITH BLACK ABSORBER CONDITIONS
   GEODST_Int , POINTER :: NZC(:)    ! NZONE              ! ZONE CLASSIFICATIONS
   GEODST_Int , POINTER :: NZNR(:)   ! NREG               ! ZONE NUMBER ASSIGNED TO EACH REGION

   ! CARD TYPE 6           ---REGION ASSIGNMENTS TO COARSE MESH INTERVALS (6D RECORD)---
   ! PRESENT IF IGOM.GT.0 AND NRASS.EQ.0
   ! CARD TYPE 7           ---REGION ASSIGNMENTS TO FINE MESH INTERVALS (7D RECORD)---
   ! PRESENT IF IGOM.GT.0 AND NRASS.EQ.1
   GEODST_Int , DIMENSION(:,:,:), POINTER :: MR  ! (6) NCINTI,NCINTJ,NCINTK ! REGION NUMBERS ASSIGNED TO COARSE MESH INTERVALS
                                                             ! (7) NINTI,NINTJ,NINTK    ! REGION NUMBERS ASSIGNED TO FINE MESH INTERVALS
END TYPE

DIF3D_Int            :: MODULE_OUT = 6                  ! Module output unit
PRIVATE MODULE_OUT
TYPE (GEODST_DATA), SAVE :: MODULE_GEODST  ! A optional use "common block" data type for the module

CONTAINS

!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_ASSIGNPRINTINFO(MODULE_OUTPUT_UNIT)
!  Sets the output unit for the module
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GEODST_ASSIGNPRINTINFO(MODULE_OUTPUT_UNIT)
DIF3D_Int MODULE_OUTPUT_UNIT    ! Reset the object output unit

MODULE_OUT = MODULE_OUTPUT_UNIT

END SUBROUTINE GEODST_ASSIGNPRINTINFO

!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_DEFINE(User_GEODST,USER_IGOM,USER_NZONE,USER_NREG,USER_NCINTI,USER_NCINTJ,USER_NCINTK,        &
!                     USER_NINTI,USER_NINTJ,USER_NINTK,USER_NBS,USER_NBCS,USER_NIBCS,USER_NZWBB,USER_NRASS   )
!  Defines a GEODST_DATA type based upon user input
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GEODST_DEFINE(User_GEODST,USER_IGOM,USER_NZONE,USER_NREG,USER_NCINTI,USER_NCINTJ,USER_NCINTK,        &
                              USER_NINTI,USER_NINTJ,USER_NINTK,USER_NBS,USER_NBCS,USER_NIBCS,USER_NZWBB,USER_NRASS   )
IMPLICIT NONE
! Passed in
TYPE (GEODST_DATA) User_GEODST                  ! A user GEODST data structure
DIF3D_Int USER_IGOM       ! GEOMETRY TYPE 
DIF3D_Int USER_NZONE      ! NUMBER OF ZONES (EACH HOMOGENEOUS IN NEUTRONICS PROBLEM - A ZONE CONTAINS ONE OR MORE REGIONS)
DIF3D_Int USER_NREG       ! NUMBER OF REGIONS
DIF3D_Int USER_NCINTI     ! NUMBER OF FIRST  DIMENSION COARSE MESH INTERVALS
DIF3D_Int USER_NCINTJ     ! NUMBER OF SECOND DIMENSION COARSE MESH INTERVALS NCINTJ = 1 FOR ONE DIMENSIONAL CASE
DIF3D_Int USER_NCINTK     ! NUMBER OF THIRD  DIMENSION COARSE MESH INTERVALS NCINTK = 1 FOR ONE AND TWO DIMENSIONAL CASES
DIF3D_Int USER_NINTI      ! NUMBER OF FIRST  DIMENSION FINE   MESH INTERVALS
DIF3D_Int USER_NINTJ      ! NUMBER OF SECOND DIMENSION FINE   MESH INTERVALS NINTJ = 1  FOR ONE DIMENSIONAL CASE
DIF3D_Int USER_NINTK      ! NUMBER OF THIRD  DIMENSION FINE   MESH INTERVALS NINTK = 1  FOR ONE AND TWO DIMENSION CASES
DIF3D_Int USER_NBS        ! NUMBER OF BUCKLING SPECIFICATIONS
DIF3D_Int USER_NBCS       ! NUMBER OF CONSTANTS FOR EXTERNAL BOUNDARIES
DIF3D_Int USER_NIBCS      ! NUMBER OF CONSTANTS FOR INTERNAL BOUNDARIES
DIF3D_Int USER_NZWBB      ! NUMBER OF ZONES WHICH ARE BLACK ABSORBERS
DIF3D_Int USER_NRASS      ! REGION ASSIGNMENTS
! Local junk
DIF3D_Int IOS,I

100 FORMAT('[GEODST]...SORRY, BUT I MUST STOP')
101 FORMAT('[GEODST]',107('.'))
102 FORMAT('')
105 FORMAT('[GEODST]...There was a fatal error that occured in (DEFINE)',56('.'))

IF (User_GEODST%DEFINED) CALL GEODST_VOID(User_GEODST)   ! If the variable is already defined, void it

#ifdef DIF3D_Debug
   WRITE(MODULE_OUT,'("[GEODST]...GEOMETRY TYPE...................................",40("."),I16)') USER_IGOM
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF ZONES.................................",40("."),I16)') USER_NZONE
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF REGIONS...............................",40("."),I16)') USER_NREG
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF FIRST  DIMENSION COARSE MESH INTERVALS",40("."),I16)') USER_NCINTI
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF SECOND DIMENSION COARSE MESH INTERVALS",40("."),I16)') USER_NCINTJ
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF THIRD  DIMENSION COARSE MESH INTERVALS",40("."),I16)') USER_NCINTK
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF FIRST  DIMENSION FINE   MESH INTERVALS",40("."),I16)') USER_NINTI
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF SECOND DIMENSION FINE   MESH INTERVALS",40("."),I16)') USER_NINTJ
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF THIRD  DIMENSION FINE   MESH INTERVALS",40("."),I16)') USER_NINTK
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF BUCKLING SPECIFICATIONS...............",40("."),I16)') USER_NBS
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF CONSTANTS FOR EXTERNAL BOUNDARIES.....",40("."),I16)') USER_NBCS
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF CONSTANTS FOR INTERNAL BOUNDARIES.....",40("."),I16)') USER_NIBCS
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF ZONES WHICH ARE BLACK ABSORBERS.......",40("."),I16)') USER_NZWBB
   WRITE(MODULE_OUT,'("[GEODST]...REGION ASSIGNMENT...............................",40("."),I16)') USER_NRASS
#endif

IF ((USER_IGOM .LT. 0) .OR. (USER_IGOM .GT. 18)) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC GEOMETRY TYPE",63("."),I16)') USER_IGOM
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (USER_NZONE .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC NUMBER OF ZONES",61("."),I16)') USER_NZONE
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (USER_NREG .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC NUMBER OF REGIONS",59("."),I16)') USER_NREG
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (USER_NCINTI .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC NUMBER OF FIRST DIMENSION COARSE MESH INTERVALS",29("."),I16)') USER_NCINTI
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (USER_NCINTJ .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC NUMBER OF SECOND DIMENSION COARSE MESH INTERVALS",28("."),I16)') USER_NCINTJ
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (USER_NCINTK .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC NUMBER OF THIRD DIMENSION COARSE MESH INTERVALS",29("."),I16)') USER_NCINTK
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (USER_NINTI .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC NUMBER OF FIRST DIMENSION FINE MESH INTERVALS",31("."),I16)') USER_NINTI
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (USER_NINTJ .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC NUMBER OF SECOND DIMENSION FINE MESH INTERVALS",30("."),I16)') USER_NINTJ
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (USER_NINTK .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC NUMBER OF THIRD DIMENSION FINE MESH INTERVALS",31("."),I16)') USER_NINTK
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (USER_NBS .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC NUMBER OF BUCKLING SPECIFICATIONS",43("."),I16)') USER_NBS
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (USER_NBCS .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC NUMBER OF EXTERNAL BOUNDARY CONSTANTS",39("."),I16)') USER_NBCS
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (USER_NIBCS .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC NUMBER OF INTERNAL BOUNDARY CONSTANTS",39("."),I16)') USER_NIBCS
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF ((USER_NZWBB .LT. 0) .OR. (USER_NZWBB .GT. USER_NZONE)) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC NUMBER OF BLACK ABSORBERS ZONES",45("."),I16)') USER_NZWBB
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF ((USER_NRASS .NE. 0) .AND. (USER_NRASS .NE. 1)) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...UNREALISTIC REGION ASSIGNMENT",59("."),I16)') USER_NRASS
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

! Allocate component arrays
ALLOCATE(User_GEODST%XMESH(USER_NCINTI+1),User_GEODST%IFINTS(USER_NCINTI),                               &
         User_GEODST%YMESH(USER_NCINTJ+1),User_GEODST%JFINTS(USER_NCINTJ),                               &
         User_GEODST%ZMESH(USER_NCINTK+1),User_GEODST%KFINTS(USER_NCINTK),                               &
         User_GEODST%VOLR(USER_NREG),User_GEODST%BSQ(USER_NBS),User_GEODST%BNDC(USER_NBCS),User_GEODST%BNCI(USER_NIBCS),   &
         User_GEODST%NZHBB(USER_NZWBB),User_GEODST%NZC(USER_NZONE),User_GEODST%NZNR(USER_NREG),                   &
         STAT = IOS)
IF (IOS .NE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...INSUFFICIENT MEMORY FOR PRIMARIES OR PREVIOUSLY DEFINED ARRAYS",42("."))')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

! Allocate the secondary component array
IF ((USER_IGOM .GT. 0) .AND. (USER_NRASS .EQ. 0)) ALLOCATE(User_GEODST%MR(USER_NCINTI,USER_NCINTJ,USER_NCINTK),STAT = IOS)
IF ((USER_IGOM .GT. 0) .AND. (USER_NRASS .EQ. 1)) ALLOCATE(User_GEODST%MR(USER_NINTI,USER_NINTJ,USER_NINTK),STAT = IOS)
IF (IOS .NE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...INSUFFICIENT MEMORY FOR MR ARRAY OR PREVIOUSLY DEFINED ARRAY",44("."))')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

! Initialize the arrays
User_GEODST%XMESH  = 0.0d0
User_GEODST%IFINTS = 0
User_GEODST%YMESH  = 0.0d0
User_GEODST%JFINTS = 0
User_GEODST%ZMESH  = 0.0d0
User_GEODST%KFINTS = 0
User_GEODST%VOLR   = 0.0d0
IF (USER_NBS .GT. 0)   User_GEODST%BSQ  = 0.0d0
IF (USER_NBCS .GT. 0)  User_GEODST%BNDC = 0.0d0
IF (USER_NIBCS .GT. 0) User_GEODST%BNCI = 0.0d0
IF (USER_NZWBB .GT. 0) User_GEODST%NZHBB = 0
User_GEODST%NZC    = 0
User_GEODST%NZNR   = 0
IF (USER_IGOM .GT. 0) User_GEODST%MR = 0

! Finish
User_GEODST%DEFINED = .TRUE.
User_GEODST%IGOM   = USER_IGOM
User_GEODST%NZONE  = USER_NZONE
User_GEODST%NREG   = USER_NREG
User_GEODST%NCINTI = USER_NCINTI
User_GEODST%NCINTJ = USER_NCINTJ
User_GEODST%NCINTK = USER_NCINTK
User_GEODST%NINTI  = USER_NINTI
User_GEODST%NINTJ  = USER_NINTJ
User_GEODST%NINTK  = USER_NINTK
User_GEODST%NBS    = USER_NBS
User_GEODST%NBCS   = USER_NBCS
User_GEODST%NIBCS  = USER_NIBCS
User_GEODST%NZWBB  = USER_NZWBB
User_GEODST%NRASS  = USER_NRASS

END SUBROUTINE GEODST_DEFINE

!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_Copy(Source_GEODST,Destination_GEODST)
!  Makes a copy of an existing GEODST data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GEODST_Copy(Source_GEODST,Destination_GEODST)
IMPLICIT NONE
! Passed in
TYPE (GEODST_DATA) Source_GEODST,Destination_GEODST                  ! A user GEODST data structure
! Local variables
DIF3D_Int IGOM,NZONE,NREG,NCINTI,NCINTJ,NCINTK,NINTI,NINTJ,NINTK,NBS,NBCS,NIBCS,NZWBB,NRASS
! Local junk
DIF3D_Int IOS,I

IF (Destination_GEODST%DEFINED) CALL GEODST_VOID(Destination_GEODST)   ! If the variable is already defined, void it

IF (Source_GEODST%Defined) THEN
   IGOM   = Source_GEODST%IGOM  
   NZONE  = Source_GEODST%NZONE 
   NREG   = Source_GEODST%NREG  
   NCINTI = Source_GEODST%NCINTI
   NCINTJ = Source_GEODST%NCINTJ
   NCINTK = Source_GEODST%NCINTK
   NINTI  = Source_GEODST%NINTI 
   NINTJ  = Source_GEODST%NINTJ 
   NINTK  = Source_GEODST%NINTK 
   NBS    = Source_GEODST%NBS   
   NBCS   = Source_GEODST%NBCS  
   NIBCS  = Source_GEODST%NIBCS 
   NZWBB  = Source_GEODST%NZWBB 
   NRASS  = Source_GEODST%NRASS 
   CALL GEODST_DEFINE(Destination_GEODST,IGOM,NZONE,NREG,NCINTI,NCINTJ,NCINTK,NINTI,NINTJ,NINTK,NBS,NBCS,NIBCS,NZWBB,NRASS)
   Destination_GEODST%HNAME   = Source_GEODST%HNAME
   Destination_GEODST%HUSE    = Source_GEODST%HUSE
   Destination_GEODST%NZCL    = Source_GEODST%NZCL
   Destination_GEODST%IMB1    = Source_GEODST%IMB1
   Destination_GEODST%IMB2    = Source_GEODST%IMB2
   Destination_GEODST%JMB1    = Source_GEODST%JMB1
   Destination_GEODST%JMB2    = Source_GEODST%JMB2
   Destination_GEODST%KMB1    = Source_GEODST%KMB1
   Destination_GEODST%KMB2    = Source_GEODST%KMB2
   Destination_GEODST%NTRIAG  = Source_GEODST%NTRIAG
   Destination_GEODST%NTHPT   = Source_GEODST%NTHPT
   Destination_GEODST%NGOP    = Source_GEODST%NGOP
   Destination_GEODST%XMESH   = Source_GEODST%XMESH
   Destination_GEODST%IFINTS  = Source_GEODST%IFINTS
   Destination_GEODST%YMESH   = Source_GEODST%YMESH
   Destination_GEODST%JFINTS  = Source_GEODST%JFINTS
   Destination_GEODST%ZMESH   = Source_GEODST%ZMESH
   Destination_GEODST%KFINTS  = Source_GEODST%KFINTS
   Destination_GEODST%VOLR    = Source_GEODST%VOLR
   Destination_GEODST%BSQ     = Source_GEODST%BSQ 
   Destination_GEODST%BNDC    = Source_GEODST%BNDC
   Destination_GEODST%BNCI    = Source_GEODST%BNCI
   Destination_GEODST%NZHBB   = Source_GEODST%NZHBB
   Destination_GEODST%NZC     = Source_GEODST%NZC
   Destination_GEODST%NZNR    = Source_GEODST%NZNR
   Destination_GEODST%MR      = Source_GEODST%MR
END IF

END SUBROUTINE GEODST_Copy

!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_NRASS1(User_GEODST)
!  Converts the coarse mesh mapping (NRASS=0) to fine mesh mapping (NRASS1) in GEODST_DATA
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GEODST_NRASS1(User_GEODST)
IMPLICIT NONE
! Passed in
TYPE (GEODST_DATA) User_GEODST                  ! A user GEODST data structure
! Local variables
TYPE (GEODST_DATA) Temp_GEODST                 ! A user GEODST data structure
! Local variables
DIF3D_Int IGOM,NZONE,NREG,NCINTI,NCINTJ,NCINTK,NINTI,NINTJ,NINTK,NBS,NBCS,NIBCS,NZWBB,NRASS
! Local junk
DIF3D_Int IOS,I,J,K,II,JJ,KK,III,JJJ,KKK,IK,IJ
DIF3D_Real TempReal

100 FORMAT('[GEODST]...SORRY, BUT I MUST STOP')
101 FORMAT('[GEODST]',107('.'))
102 FORMAT('')
105 FORMAT('[GEODST]...There was a fatal error that occured in (DEFINE)',56('.'))

IF ((User_GEODST%Defined) .AND. (User_GEODST%NRASS .EQ. 0)) THEN
   !WRITE(MODULE_OUT,*)'Data before altering'
   !CALL GEODST_Print(User_GEODST)
   IGOM   = User_GEODST%IGOM  
   NZONE  = User_GEODST%NZONE 
   NREG   = User_GEODST%NREG  
   NCINTI = User_GEODST%NCINTI
   NCINTJ = User_GEODST%NCINTJ
   NCINTK = User_GEODST%NCINTK
   NINTI  = User_GEODST%NINTI 
   NINTJ  = User_GEODST%NINTJ 
   NINTK  = User_GEODST%NINTK 
   NBS    = User_GEODST%NBS   
   NBCS   = User_GEODST%NBCS  
   NIBCS  = User_GEODST%NIBCS 
   NZWBB  = User_GEODST%NZWBB 
   NRASS  = 1 ! As opposed to 0
   CALL GEODST_DEFINE(Temp_GEODST,IGOM,NZONE,NREG,NINTI,NINTJ,NINTK,NINTI,NINTJ,NINTK,NBS,NBCS,NIBCS,NZWBB,NRASS)
   ! Copy the rest of the data
   Temp_GEODST%HNAME   = User_GEODST%HNAME
   Temp_GEODST%HUSE    = User_GEODST%HUSE
   Temp_GEODST%NZCL    = User_GEODST%NZCL
   Temp_GEODST%IMB1    = User_GEODST%IMB1
   Temp_GEODST%IMB2    = User_GEODST%IMB2
   Temp_GEODST%JMB1    = User_GEODST%JMB1
   Temp_GEODST%JMB2    = User_GEODST%JMB2
   Temp_GEODST%KMB1    = User_GEODST%KMB1
   Temp_GEODST%KMB2    = User_GEODST%KMB2
   Temp_GEODST%NTRIAG  = User_GEODST%NTRIAG
   Temp_GEODST%NTHPT   = User_GEODST%NTHPT
   Temp_GEODST%NGOP    = User_GEODST%NGOP
                         
   Temp_GEODST%VOLR(1:NREG)   = User_GEODST%VOLR(1:NREG)
   Temp_GEODST%BSQ(1:NBS)     = User_GEODST%BSQ(1:NBS)
   Temp_GEODST%BNDC(1:NBCS)   = User_GEODST%BNDC(1:NBCS)
   Temp_GEODST%BNCI(1:NIBCS)  = User_GEODST%BNCI(1:NIBCS)
   Temp_GEODST%NZHBB(1:NZWBB) = User_GEODST%NZHBB(1:NZWBB)
   Temp_GEODST%NZC(1:NZONE)   = User_GEODST%NZC(1:NZONE)
   Temp_GEODST%NZNR(1:NREG)   = User_GEODST%NZNR(1:NREG)
   ! Convert the meshing details
   IF ((User_GEODST%IGOM .GE.  1) .AND. (User_GEODST%IGOM .LE.  2)  .OR. &
       (User_GEODST%IGOM .GE.  6) .AND. (User_GEODST%IGOM .LE.  8)  .OR. &
       (User_GEODST%IGOM .GE. 11) .AND. (User_GEODST%IGOM .LE. 16)       )  THEN
      III = 1
      DO I = 1,NCINTI
         TempReal = User_GEODST%IFINTS(I)
         TempReal = (User_GEODST%XMESH(I+1)-User_GEODST%XMESH(I))/TempReal
         DO II = 1,User_GEODST%IFINTS(I)
            III = III + 1
            Temp_GEODST%XMESH(III) = User_GEODST%XMESH(I) + TempReal * II
            Temp_GEODST%IFINTS(III-1) = 1
         END DO
      END DO
      JJJ = 1
      DO J = 1,NCINTJ
         TempReal = User_GEODST%JFINTS(J)
         TempReal = (User_GEODST%YMESH(J+1)-User_GEODST%YMESH(J))/TempReal
         DO JJ = 1,User_GEODST%JFINTS(J)
            JJJ = JJJ + 1
            Temp_GEODST%YMESH(JJJ) = User_GEODST%YMESH(J) + TempReal * JJ
            Temp_GEODST%JFINTS(JJJ-1) = 1
         END DO
      END DO
   ELSE
      Temp_GEODST%XMESH = User_GEODST%XMESH
      Temp_GEODST%YMESH = User_GEODST%YMESH
      Temp_GEODST%IFINTS = 1
      Temp_GEODST%JFINTS = 1
   END IF
   KKK = 1
   DO K = 1,NCINTK
      TempReal = User_GEODST%KFINTS(K)
      TempReal = (User_GEODST%ZMESH(K+1)-User_GEODST%ZMESH(K))/TempReal
      DO KK = 1,User_GEODST%KFINTS(K)
         KKK = KKK + 1
         Temp_GEODST%ZMESH(KKK) = User_GEODST%ZMESH(K) + TempReal * KK
         Temp_GEODST%KFINTS(KKK-1) = 1
      END DO
   END DO
   ! Convert MR(NCINTI,NCINTJ,NCINTK) to MR(NINTI,NINTJ,NINTK)
   KKK = 0
   DO K = 1,NCINTK
      IK = User_GEODST%KFINTS(K)
      IF (IK .EQ. 0) IK = 1
      DO KK = 1,IK
         KKK = KKK + 1
         JJJ = 0
         DO J = 1,NCINTJ
            IJ = User_GEODST%JFINTS(J)
            IF (IJ .EQ. 0) IJ = 1
            DO JJ = 1,IJ
               JJJ = JJJ + 1
               III = 0
               DO I = 1,NCINTI
                  DO II = 1,User_GEODST%IFINTS(I)
                     III = III + 1
                     Temp_GEODST%MR(III,JJJ,KKK) = User_GEODST%MR(I,J,K)
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
   ! Copy the data back
   CALL GEODST_Copy(Temp_GEODST,User_GEODST)
   CALL GEODST_VOID(Temp_GEODST)
   !WRITE(MODULE_OUT,*)'Data after altering'
   !CALL GEODST_Print(User_GEODST)
END IF

END SUBROUTINE GEODST_NRASS1

!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_AddMeshes(User_GEODST,IJK,NumBoundaries,MeshBoundaries)
!  This routine adds the user specified meshes to an existing domain
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GEODST_AddMeshes(User_GEODST,IJK,NumBoundaries,MeshBoundaries)
IMPLICIT NONE
! Passed in
TYPE (GEODST_DATA) User_GEODST                  ! A user GEODST data structure
DIF3D_Int  IJK    ! pick 1,2,3 for I,J,K
DIF3D_Int  NumBoundaries
DIF3D_Real MeshBoundaries(NumBoundaries)
! Local variables
TYPE (GEODST_DATA) Temp_GEODST                 ! A user GEODST data structure
! Local variables
DIF3D_Int IGOM,NZONE,NREG,NCINTI,NCINTJ,NCINTK,NINTI,NINTJ,NINTK,NBS,NBCS,NIBCS,NZWBB,NRASS
! Local junk
DIF3D_Int IOS,I,J,K,II,JJ,KK,III,JJJ,KKK
DIF3D_Real TempReal
DIF3D_Int  IBounds
DIF3D_Int  BoundIndex(10000)
DIF3D_Real Boundaries(10000)

100 FORMAT('[GEODST]...SORRY, BUT I MUST STOP')
101 FORMAT('[GEODST]',107('.'))
102 FORMAT('')
105 FORMAT('[GEODST]...There was a fatal error that occured in (AddMeshes)',56('.'))

IF (NumBoundaries .EQ. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...Invalid number of axial boundaries ",I16)') NumBoundaries
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF ((IJK .LT. 1) .OR. (IJK .GT. 3)) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...Invalid ijk selection value ",I16)') IJK
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

! Check IJK
I = 0
IF (User_GEODST%IGOM .LT. 6) THEN ! 1D
   IF (IJK .GT. 1) I = 1
ELSE IF (User_GEODST%IGOM .LT. 12) THEN ! 2D
   IF (IJK .GT. 2) I = 1
END IF
IF (I .EQ. 1) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...Geometry + IJK selection is not possible ",I16)') IJK
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (.NOT. User_GEODST%Defined) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...Structure was not defined?")')
   WRITE(MODULE_OUT,100)
   CALL Abort
ELSE
   ! Convert the input to NRASS=1
   IF (User_GEODST%NRASS .EQ. 0) CALL GEODST_NRASS1(User_GEODST)
   IGOM   = User_GEODST%IGOM  
   NZONE  = User_GEODST%NZONE 
   NREG   = User_GEODST%NREG  
   NCINTI = User_GEODST%NINTI
   NCINTJ = User_GEODST%NINTJ
   NCINTK = User_GEODST%NINTK
   NINTI  = User_GEODST%NINTI 
   NINTJ  = User_GEODST%NINTJ 
   NINTK  = User_GEODST%NINTK 
   NBS    = User_GEODST%NBS   
   NBCS   = User_GEODST%NBCS  
   NIBCS  = User_GEODST%NIBCS 
   NZWBB  = User_GEODST%NZWBB 
   NRASS  = 1 ! As opposed to 0
   ! Copy the existing data
   IF (IJK .EQ. 1) THEN
      III = NCINTI
      IBounds = NCINTI ! The new number of meshes to add
      DO I = 1,IBounds+1
         BoundIndex(I) = I
         Boundaries(I) = User_GEODST%XMESH(I)
      END DO
   ELSE IF (IJK .EQ. 2) THEN
      III = NCINTJ
      IBounds = NCINTJ ! The new number of meshes to add
      DO I = 1,IBounds+1
         BoundIndex(I) = I
         Boundaries(I) = User_GEODST%YMESH(I)
      END DO
   ELSE
      III = NCINTK
      IBounds = NCINTK ! The new number of meshes to add
      DO I = 1,IBounds+1
         BoundIndex(I) = I
         Boundaries(I) = User_GEODST%ZMESH(I)
      END DO
   END IF
   ! Add the new boundaries in
   DO II = 1,NumBoundaries
      I = 0
      DO I = 1,IBounds
         IF (Boundaries(I) .EQ. MeshBoundaries(II)) THEN ! Nothing to do
            EXIT
         ELSE IF ( (MeshBoundaries(II) .GT. Boundaries(I)  ) .AND. &
                   (MeshBoundaries(II) .LT. Boundaries(I+1))       ) THEN
            ! Move the existing data upwards
            DO J = IBounds+1,I+1,-1
               BoundIndex(J+1) = BoundIndex(J)
               Boundaries(J+1) = Boundaries(J)
            END DO
            IBounds = IBounds + 1 ! We add the new mesh point between I and I+1
            BoundIndex(I+1) = BoundIndex(I) ! Duplicate what was at position I
            Boundaries(I+1) = MeshBoundaries(II)
            EXIT ! Nothing more to do
         END IF
      END DO
   END DO
   ! Resize the data
   IF (IBounds .NE. III) THEN
      IF      (IJK .EQ. 1) THEN
         CALL GEODST_DEFINE(Temp_GEODST,IGOM,NZONE,NREG,IBounds,NCINTJ,NCINTK,IBounds,NINTJ,NINTK,NBS,NBCS,NIBCS,NZWBB,NRASS)
      ELSE IF (IJK .EQ. 2) THEN
         CALL GEODST_DEFINE(Temp_GEODST,IGOM,NZONE,NREG,NCINTI,IBounds,NCINTK,NINTI,IBounds,NINTK,NBS,NBCS,NIBCS,NZWBB,NRASS)
      ELSE
         CALL GEODST_DEFINE(Temp_GEODST,IGOM,NZONE,NREG,NCINTI,NCINTJ,IBounds,NINTI,NINTJ,IBounds,NBS,NBCS,NIBCS,NZWBB,NRASS)
      END IF
   END IF
   ! Update the data structure
   IF (Temp_GEODST%Defined) THEN
      Temp_GEODST%FILENAME = User_GEODST%FILENAME
      Temp_GEODST%HNAME    = User_GEODST%HNAME
      Temp_GEODST%HUSE     = User_GEODST%HUSE
      Temp_GEODST%IVERS    = User_GEODST%IVERS 
      Temp_GEODST%NZCL     = User_GEODST%NZCL   
      Temp_GEODST%IMB1     = User_GEODST%IMB1   
      Temp_GEODST%IMB2     = User_GEODST%IMB2   
      Temp_GEODST%JMB1     = User_GEODST%JMB1   
      Temp_GEODST%JMB2     = User_GEODST%JMB2   
      Temp_GEODST%KMB1     = User_GEODST%KMB1   
      Temp_GEODST%KMB2     = User_GEODST%KMB2   
      Temp_GEODST%NTRIAG   = User_GEODST%NTRIAG 
      Temp_GEODST%NTHPT    = User_GEODST%NTHPT  
      Temp_GEODST%NGOP     = User_GEODST%NGOP  
      Temp_GEODST%VOLR     = User_GEODST%VOLR
      Temp_GEODST%BSQ      = User_GEODST%BSQ
      Temp_GEODST%BNDC     = User_GEODST%BNDC
      Temp_GEODST%BNCI     = User_GEODST%BNCI
      Temp_GEODST%NZHBB    = User_GEODST%NZHBB
      Temp_GEODST%NZC      = User_GEODST%NZC
      Temp_GEODST%NZNR     = User_GEODST%NZNR
      Temp_GEODST%IFINTS   = 1
      Temp_GEODST%JFINTS   = 1
      Temp_GEODST%KFINTS   = 1
      IF (IJK .EQ. 1) THEN
         Temp_GEODST%XMESH(1:IBounds+1) = Boundaries(1:IBounds+1)
         Temp_GEODST%YMESH              = User_GEODST%YMESH
         Temp_GEODST%ZMESH              = User_GEODST%ZMESH
         DO K = 1,NCINTK
            DO J = 1,NCINTJ
               DO I = 1,IBounds
                  II = BoundIndex(I)
                  Temp_GEODST%MR(I,J,K) = User_GEODST%MR(II,J,K)
               END DO
            END DO
         END DO
      ELSE IF (IJK .EQ. 2) THEN
         Temp_GEODST%XMESH              = User_GEODST%XMESH
         Temp_GEODST%YMESH(1:IBounds+1) = Boundaries(1:IBounds+1)
         Temp_GEODST%ZMESH              = User_GEODST%ZMESH
         DO K = 1,NCINTK
            DO I = 1,IBounds
               II = BoundIndex(I)
               DO J = 1,NCINTI
                  Temp_GEODST%MR(J,I,K) = User_GEODST%MR(J,II,K)
               END DO
            END DO
         END DO
      ELSE
         Temp_GEODST%XMESH              = User_GEODST%XMESH
         Temp_GEODST%YMESH              = User_GEODST%YMESH
         Temp_GEODST%ZMESH(1:IBounds+1) = Boundaries(1:IBounds+1)
               DO I = 1,IBounds
                  II = BoundIndex(I)
         DO J = 1,NCINTJ
            DO K = 1,NCINTI
                  Temp_GEODST%MR(K,J,I) = User_GEODST%MR(K,J,II)
            END DO
         END DO
               END DO
      END IF
      CALL GEODST_Copy(Temp_GEODST,User_GEODST) ! Overwrite the old data
      CALL GEODST_Void(Temp_GEODST)
   END IF
END IF

END SUBROUTINE GEODST_AddMeshes

!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_LocatePoint(User_GEODST,X,Y,Z,oINode,oJNode,oKNode,o_Reg,oX,oY,oZ)
!  This routine identifies the GEODST I,J,K position and local coordinates of a given point within that mesh
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GEODST_LocatePoint(User_GEODST,X,Y,Z,oINode,oJNode,oKNode,o_Reg,oX,oY,oZ)
IMPLICIT NONE
! Passed in
TYPE (GEODST_DATA) User_GEODST                  ! A user GEODST data structure
DIF3D_Real X,Y,Z    ! pick 1,2,3 for I,J,K
DIF3D_Int  oINode,oJNode,oKNode,o_Reg
DIF3D_Real oX,oY,oZ

! Local junk
DIF3D_Int IOS,I,J,K,II,JJ,KK,III,JJJ,KKK,HexOption
DIF3D_Real TempReal
DIF3D_Real ::  cos30=0.866025403784439d0
DIF3D_Real :: icos30=1.154700538379251d0
DIF3D_Real :: SomeReal,Delta,Rmin,HexPitch,DiagPitch,HalfPitch,GeomScaling
DIF3D_Real :: xcenter,ycenter,XX,YY
DIF3D_Log  :: HitSomething

100 FORMAT('[GEODST]...SORRY, BUT I MUST STOP')
101 FORMAT('[GEODST]',107('.'))
102 FORMAT('')
105 FORMAT('[GEODST]...There was a fatal error that occured in (LocatePoint)',56('.'))

IF (.NOT. User_GEODST%Defined) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...Data structure not defined ")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

oINode = 0
oJNode = 0
oKNode = 0
o_Reg  = 0
oX = 0.0d0
oY = 0.0d0
oZ = 0.0d0

! GEOMETRY TYPE 
!       0- POINT (FUNDAMENTAL MODE),  1- SLAB  ,  2- CYLINDER  ,  3- SPHERE               ,
!       6- X-Y                     ,  7- R-Z   ,  8- THETA-R   ,  9- UNIFORM TRIANGULAR   ,
!      10- HEXAGONAL (1 MESH POINT IN EACH HEXAGONAL ELEMENT)  , 11- R-THETA              ,
!      12- R-THETA-Z,              , 13- R-THETA-ALPHA         , 14- X-Y-Z                ,
!      15- THETA-R-Z,              , 16- THETA-R-ALPHA,        , 17- UNIFORM TRIANGULAR-Z ,
!      18- HEXAGON-Z   (MESH POINTS AS IN 10 ABOVE) 

! There is a different algorithm for each geometry type
IF      ((User_GEODST%IGOM .EQ. 1) .OR. (User_GEODST%IGOM .EQ. 6) .OR. (User_GEODST%IGOM .EQ. 14)) THEN
   ! The GEODST ordering is easy to compute if we know the I,J,K mesh position
   I=0
   J=1
   K=1
   III = 0
   JJJ = 1
   KKK = 1
   HitSomething = .FALSE.
   IF ((X .GE. User_GEODST%XMESH(1)) .AND. (X .LE. User_GEODST%XMESH(User_GEODST%NCINTI+1))) THEN
      C1D: DO I = 1,User_GEODST%NCINTI
         SomeReal = User_GEODST%IFINTS(I)
         Delta = (User_GEODST%XMESH(I+1)-User_GEODST%XMESH(I))/SomeReal
         Rmin = User_GEODST%XMESH(I) - Delta
         DO II = 1,User_GEODST%IFINTS(I)
            III = III + 1
            Rmin = Rmin + Delta
            IF ((Rmin .LE. X) .AND. (X .LE. Rmin+Delta)) THEN
               HitSomething = .TRUE.
               oX = -0.5d0 + (X - Rmin)/Delta
               EXIT C1D
            END IF
         END DO
      END DO C1D
   END IF
   IF (.NOT. HitSomething) III = -1
   IF (User_GEODST%IGOM .NE. 1) THEN
      HitSomething = .FALSE.
      JJJ = 0
      IF ((Y .GE. User_GEODST%YMESH(1)) .AND. (Y .LE. User_GEODST%YMESH(User_GEODST%NCINTJ+1))) THEN
         C2D: DO J = 1,User_GEODST%NCINTJ
            SomeReal = User_GEODST%JFINTS(J)
            Delta = (User_GEODST%YMESH(J+1)-User_GEODST%YMESH(J))/SomeReal
            Rmin = User_GEODST%YMESH(J) - Delta
            DO JJ = 1,User_GEODST%JFINTS(J)
               JJJ = JJJ + 1
               Rmin = Rmin + Delta
               IF ((Rmin .LE. Y) .AND. (Y .LE. Rmin+Delta)) THEN
                  HitSomething = .TRUE.
                  oY = -0.5d0 + (Y - Rmin)/Delta
                  EXIT C2D
               END IF
            END DO
         END DO C2D
      END IF
      IF (.NOT. HitSomething) JJJ = -1
   END IF
   IF (User_GEODST%IGOM .EQ. 14) THEN
      HitSomething = .FALSE.
      KKK = 0
      IF ((Z .GE. User_GEODST%ZMESH(1)) .AND. (Z .LE. User_GEODST%ZMESH(User_GEODST%NCINTK+1))) THEN
         C3D: DO K = 1,User_GEODST%NCINTK
            SomeReal = User_GEODST%KFINTS(K)
            Delta = (User_GEODST%ZMESH(K+1)-User_GEODST%ZMESH(K))/SomeReal
            Rmin = User_GEODST%ZMESH(K) - Delta
            DO KK = 1,User_GEODST%KFINTS(K)
               KKK = KKK + 1
               Rmin = Rmin + Delta
               IF ((Rmin .LE. Z) .AND. (Z .LE. Rmin+Delta)) THEN
                  HitSomething = .TRUE.
                  oZ = -0.5d0 + (Z - Rmin)/Delta
                  EXIT C3D
               END IF
            END DO
         END DO C3D
      END IF
      IF (.NOT. HitSomething) KKK = -1
   END IF
   IF ((III .GT. 0) .AND. (JJJ .GT. 0) .AND. (KKK .GT. 0)) THEN
      oINode = III
      oJNode = JJJ
      IF (User_GEODST%NRASS .EQ. 0) THEN
         o_Reg = User_GEODST%MR(I,J,K)
      ELSE
         o_Reg = User_GEODST%MR(III,JJJ,KKK)
      END IF
   END IF
   oKNode  = KKK
ELSE IF ((User_GEODST%IGOM .EQ. 10) .OR. (User_GEODST%IGOM .EQ. 18)) THEN
   K = 1
   KKK = 1
   IF (User_GEODST%IGOM .EQ. 18) THEN
      KKK = 0
      IF ((Z .GE. User_GEODST%ZMESH(1)) .AND. (Z .LE. User_GEODST%ZMESH(User_GEODST%NCINTK+1))) THEN
         H3D: DO K = 1,User_GEODST%NCINTK
            SomeReal = User_GEODST%KFINTS(K)
            Delta = (User_GEODST%ZMESH(K+1)-User_GEODST%ZMESH(K))/SomeReal
            Rmin = User_GEODST%ZMESH(K) - Delta
            DO KK = 1,User_GEODST%KFINTS(K)
               KKK = KKK + 1
               Rmin = Rmin + Delta
               IF ((Rmin .LE. Z) .AND. (Z .LE. Rmin+Delta)) THEN
                  oZ = -0.5d0 + (Z - Rmin)/Delta
                  EXIT H3D
               END IF
            END DO
         END DO H3D
      END IF
   END IF
   ! I need to detect whether this is 
   ! HexOption = 0 is a full core
   ! HexOption = 1 is 120 symmetry/periodic which will have Base_GEODST%IMB1=4 and Base_GEODST%NTRIAG = 0
   ! HexOption = 2 is  60 symmetry/periodic which will have Base_GEODST%IMB1=4 and Base_GEODST%NTRIAG = 1
   HexOption = 0
   IF (User_GEODST%IMB1 .EQ. 4) THEN
      HexOption = 1
      IF (User_GEODST%NTRIAG .EQ. 1) HexOption = 2
   END IF
   !WRITE(MODULE_OUT,*)'HexOption = ',HexOption
   HexPitch  = User_GEODST%XMESH(2)
   HalfPitch = 0.5d0 * HexPitch
   DiagPitch = icos30 * HalfPitch
   ! The pitch of a volume=1 assembly is 1.074569931823542
   ! I need to scale the existing assembly by 1.074569931823542/Pitch in the X and Y direction
   GeomScaling = 1.074569931823542d0 / HexPitch
   IF (HexOption .EQ. 0) THEN
      ! The geometry center is at hex NINTI/2+1,NINTJ/2+1 where J+ movement is X- movement
      III = User_GEODST%NINTI/2+1
      JJJ = User_GEODST%NINTJ/2+1
      H360: DO J = 1,User_GEODST%NINTJ
         DO I = 1,User_GEODST%NINTI
            xcenter = (I-III)*HexPitch - (J-JJJ)*HalfPitch
            ycenter = (J-JJJ)*HexPitch*cos30
            !WRITE(MODULE_OUT,*)'I= ',I,' J=',J
            !WRITE(MODULE_OUT,*)'xc= ',xcenter,' yc= ', ycenter
            !WRITE(MODULE_OUT,*)'X= ',X,' Y= ',Y
            ! This transformation puts our point in the positive direction relative to the center of this hex
            XX = dabs(X - xcenter)
            YY = dabs(Y - ycenter)
            IF ((XX .GT. HalfPitch) .OR. (YY .GT. DiagPitch)) CYCLE ! Outside of the simple bounding box on the first quadrant
            ! This checks whether our point lies above or below the 60-90 degree line (should be below)
            IF (DiagPitch*HalfPitch - 0.5d0*DiagPitch*XX - HalfPitch * YY .GE. 0.0d0) THEN ! We are inside of the hex
               oX = (X - xcenter) * GeomScaling
               oY = (Y - ycenter) * GeomScaling
               EXIT H360
            END IF
         END DO
      END DO H360
   ELSE IF (HexOption .EQ. 1) THEN
      ! The center of the 1,1 hex is the geometry center where J+ movement is X- movement 
      H120: DO J = 1,User_GEODST%NINTJ
         DO I = 1,User_GEODST%NINTI
            xcenter = (I-1)*HexPitch - (J-1)*HalfPitch
            ycenter = (J-1)*HexPitch*cos30
            ! This transformation puts our point in the positive direction relative to the center of this hex
            XX = dabs(X - xcenter)
            YY = dabs(Y - ycenter)
            IF ((XX .GT. HalfPitch) .OR. (YY .GT. DiagPitch)) CYCLE ! Outside of the simple bounding box on the first quadrant
            ! This checks whether our point lies above or below the 60-90 degree line (should be below)
            IF (DiagPitch*HalfPitch - 0.5d0*DiagPitch*XX - HalfPitch * YY .GE. 0.0d0) THEN ! We are inside of the hex
               oX = (X - xcenter) * GeomScaling
               oY = (Y - ycenter) * GeomScaling
               EXIT H120
            END IF
         END DO
      END DO H120
   ELSE
      ! The center of the 1,1 hex is the geometry center where J+ movement is X+ movement 
      H60: DO J = 1,User_GEODST%NINTJ
         DO I = 1,User_GEODST%NINTI
            xcenter = (I-1)*HexPitch + (J-1)*HalfPitch
            ycenter = (J-1)*HexPitch*cos30
            !WRITE(MODULE_OUT,*)'I= ',I,' J=',J
            !WRITE(MODULE_OUT,*)'xc= ',xcenter,' yc= ', ycenter
            !WRITE(MODULE_OUT,*)'X= ',X,' Y= ',Y
            ! This transformation puts our point in the positive direction relative to the center of this hex
            XX = dabs(X - xcenter)
            YY = dabs(Y - ycenter)
            !WRITE(MODULE_OUT,*)'XX= ',XX,' YY= ',YY
            !WRITE(MODULE_OUT,*)'cX= ',HalfPitch,' cY= ',DiagPitch
            IF ((XX .GT. HalfPitch) .OR. (YY .GT. DiagPitch)) CYCLE ! Outside of the simple bounding box on the first quadrant
            !WRITE(MODULE_OUT,*)'Survived 1'
            ! This checks whether our point lies above or below the 60-90 degree line (should be below)
            !WRITE(MODULE_OUT,*)'product =',DiagPitch*HalfPitch - 0.5d0*DiagPitch*XX - HalfPitch * YY
            IF (DiagPitch*HalfPitch - 0.5d0*DiagPitch*XX - HalfPitch * YY .GE. 0.0d0) THEN ! We are inside of the hex
               oX = (X - xcenter) * GeomScaling
               oY = (Y - ycenter) * GeomScaling
               EXIT H60
            END IF
         END DO
      END DO H60
   END IF
   IF ((I .LE. User_GEODST%NINTI) .AND. (J .LE. User_GEODST%NINTJ) .AND. (KKK .GT. 0)) THEN
      oINode = I
      oJNode = J
      IF (User_GEODST%NRASS .EQ. 0) THEN
         o_Reg = User_GEODST%MR(I,J,K)
      ELSE
         o_Reg = User_GEODST%MR(I,J,KKK)
      END IF
   END IF
   oKNode  = KKK
ELSE
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...Unsupported data type, call MAS :-( ")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

END SUBROUTINE GEODST_LocatePoint

!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_ComputeCentroids(User_GEODST,Centroids)
!  This routine computes and returns the centroids of each region in the GEODST mesh
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GEODST_ComputeCentroids(User_GEODST,Centroids)
IMPLICIT NONE
! Passed in
TYPE (GEODST_DATA) User_GEODST                  ! A user GEODST data structure
DIF3D_Real Centroids(3,User_GEODST%NREG)
! Local arrays
DIF3D_Real, POINTER :: Volumes(:)
TYPE (GEODST_DATA) Copy_GEODST                  ! A user GEODST data structure
! Local junk
DIF3D_Int I,J,K,III,IOS
DIF3D_Real X_Length,Y_Length,Z_Length,Volume
DIF3D_Real :: HexC1 = 0.5D0, HexC2 = 0.866025403784439d0, HexC3 = 0.9306048591021d0
DIF3D_Real HexPitch


100 FORMAT('[GEODST]...SORRY, BUT I MUST STOP')
101 FORMAT('[GEODST]',107('.'))
102 FORMAT('')
105 FORMAT('[GEODST]...There was a fatal error that occured in (ComputeCentroids)',56('.'))

IF (.NOT. User_GEODST%Defined) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...GEODST structure must be defined in order for this to work :-/ ")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

CALL GEODST_Copy(User_GEODST,Copy_GEODST)
CALL GEODST_NRASS1(Copy_GEODST)

ALLOCATE(Volumes(Copy_GEODST%NREG),STAT=IOS)
IF (IOS .NE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...Scratch memory failure 1")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF
Volumes = 0.0D0
Centroids = 0.0d0

DO K = 1,Copy_GEODST%NINTK
   DO J = 1,Copy_GEODST%NINTJ
      DO I = 1,Copy_GEODST%NINTI
         III = Copy_GEODST%MR(I,J,K) ! The region number assigned to this mesh
         IF (III .NE. 0) THEN
            IF ((Copy_GEODST%IGOM .EQ. 1) .OR. (Copy_GEODST%IGOM .EQ. 6) .OR. (Copy_GEODST%IGOM .EQ. 14)) THEN  ! Generic cartesian x-y-z mesh option
               X_Length = Copy_GEODST%XMESH(I+1) - Copy_GEODST%XMESH(I)
               Y_Length = Copy_GEODST%YMESH(J+1) - Copy_GEODST%YMESH(J)
               Z_Length = Copy_GEODST%ZMESH(K+1) - Copy_GEODST%ZMESH(K)
               IF (Copy_GEODST%IGOM .EQ. 1) Y_Length = 1.0d0
               IF (Copy_GEODST%IGOM .LE. 6) Z_Length = 1.0d0
               Volume = X_Length*Y_Length*Z_Length
               Centroids(1,III) = Centroids(1,III) + (Copy_GEODST%XMESH(I) + 0.5d0 * X_Length)*Volume
               Centroids(2,III) = Centroids(2,III) + (Copy_GEODST%YMESH(J) + 0.5d0 * Y_Length)*Volume
               Centroids(3,III) = Centroids(3,III) + (Copy_GEODST%ZMESH(K) + 0.5d0 * Z_Length)*Volume
               Volumes(III) = Volumes(III) + Volume
            ELSE IF ((Copy_GEODST%IGOM .EQ. 10) .OR. (Copy_GEODST%IGOM .EQ. 18)) THEN ! HexZ 
               HexPitch = Copy_GEODST%XMESH(2)
               X_Length = HexPitch * HexC3
               Y_Length = X_Length
               Z_Length = Copy_GEODST%ZMESH(K+1) - Copy_GEODST%ZMESH(K)
               IF (Copy_GEODST%IGOM .EQ. 10) Z_Length = 1.0d0
               Volume = X_Length*Y_Length*Z_Length
               Centroids(1,III) = Centroids(1,III) + ((I-1)*HexPitch - (J-1)*HexC1*HexPitch) * Volume
               Centroids(2,III) = Centroids(2,III) + (                 (J-1)*HexC2*HexPitch) * Volume
               Centroids(3,III) = Centroids(3,III) + (Copy_GEODST%ZMESH(K) + 0.5d0 * Z_Length)    * Volume
               Volumes(III) = Volumes(III) + Volume
            ELSE IF (Copy_GEODST%IGOM .EQ. 17) THEN ! TriZ
               IOS = (I+1)/2   !  1,2,3,4,5,6,.... -> 1,1,2,2,3,3,4,4,...   * 2 = 2,2,4,4,6,6,8,8,...
               ! Use these as pitch
               X_Length = 0.288675134594813d0 * HexPitch
               Y_Length = 0.5d0 * HexPitch
               Z_Length = Copy_GEODST%ZMESH(K+1) - Copy_GEODST%ZMESH(K)
               Volume = X_Length*Y_Length*Z_Length
               IF (IOS*2 .EQ. I) THEN ! I was even and thus this is an upper triangle
                  Centroids(1,III) = Centroids(1,III) + (IOS-1.0d0+0.5d0+(J-1)*0.5d0)*X_Length*Volume
                  Centroids(2,III) = Centroids(2,III) + (J  -1.0d0+0.6666666667d0   )*Y_Length*Volume
               ELSE                   ! I was odd and thus this is a lower triangle
                  Centroids(1,III) = Centroids(1,III) + (IOS-1.0d0      +(J-1)*0.5d0)*X_Length*Volume
                  Centroids(2,III) = Centroids(2,III) + (J  -1.0d0+0.3333333333d0   )*Y_Length*Volume
               END IF
               Centroids(3,III) = Centroids(3,III) + (Copy_GEODST%ZMESH(K) + 0.5d0 * Z_Length)*Volume
               Volumes(III) = Volumes(III) + Volume
            ELSE IF (Copy_GEODST%IGOM .EQ. 7) THEN ! RZ
            ! This is obviously not correct as the answer should be zero for the x coordinate always :-)
            !   X_Length = Copy_GEODST%XMESH(I+1) - Copy_GEODST%XMESH(I)
            !   Y_Length = Copy_GEODST%YMESH(J+1) - Copy_GEODST%YMESH(J)
            !   Z_Length = 1.0D0 ! Copy_GEODST%ZMESH(K+1) - Copy_GEODST%ZMESH(K)
            !   Volume = (Copy_GEODST%XMESH(I+1)*Copy_GEODST%XMESH(I+1) &
            !          -  Copy_GEODST%XMESH(I)*Copy_GEODST%XMESH(I)     )*DIF3D_Definition_PI*Y_Length*Z_Length
            !   Centroids(1,III) = Centroids(1,III) + (Copy_GEODST%XMESH(I) + 0.5d0 * X_Length)*Volume
            !   Centroids(2,III) = Centroids(2,III) + (Copy_GEODST%YMESH(J) + 0.5d0 * Y_Length)*Volume
            !   Centroids(3,III) = Centroids(3,III) + (Copy_GEODST%ZMESH(K) + 0.5d0 * Z_Length)*Volume
            !   Volumes(III) = Volumes(III) + Volume
            ! For my purpose I just need this
               X_Length = Copy_GEODST%XMESH(I+1) - Copy_GEODST%XMESH(I)
               Y_Length = Copy_GEODST%YMESH(J+1) - Copy_GEODST%YMESH(J)
               Z_Length = 1.0d0
               Volume = X_Length*Y_Length*Z_Length
               Centroids(1,III) = Centroids(1,III) + (Copy_GEODST%XMESH(I) + 0.5d0 * X_Length)*Volume
               Centroids(2,III) = Centroids(2,III) + (Copy_GEODST%YMESH(J) + 0.5d0 * Y_Length)*Volume
               Centroids(3,III) = Centroids(3,III) + (Copy_GEODST%ZMESH(K) + 0.5d0 * Z_Length)*Volume
               Volumes(III) = Volumes(III) + Volume
            ELSE
               WRITE(MODULE_OUT,105)
               WRITE(MODULE_OUT,'("[GEODST]...Unknown geometry type, Call MAS :-( ")')
               WRITE(MODULE_OUT,100)
               CALL Abort
            END IF
         END IF
      END DO
   END DO
END DO
! To finish the centroid, we must divide by the volume
DO I = 1,Copy_GEODST%NREG
   Volume = Volumes(I)
   IF (Volume .NE. 0.0d0) Volume = 1.0d0/Volume
   Centroids(1,I) = Centroids(1,I) * Volume
   Centroids(2,I) = Centroids(2,I) * Volume
   Centroids(3,I) = Centroids(3,I) * Volume
END DO

DEALLOCATE(Volumes,STAT=IOS)
IF (IOS .NE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...Scratch memory failure 2")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

CALL GEODST_Void(Copy_GEODST)

END SUBROUTINE GEODST_ComputeCentroids

!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_ComputeAverageThickness(User_GEODST,Thickness)
!  This routine computes and returns the average thickness of each region in the GEODST mesh with respect to the coordinate direction
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GEODST_ComputeAverageThickness(User_GEODST,Thickness)
IMPLICIT NONE
! Passed in
TYPE (GEODST_DATA) User_GEODST                  ! A user GEODST data structure
DIF3D_Real Thickness(3,User_GEODST%NREG)
! Local arrays
DIF3D_Real, POINTER :: Volumes(:)
TYPE (GEODST_DATA) Copy_GEODST                  ! A user GEODST data structure
! Local junk
DIF3D_Int I,J,K,III,IOS
DIF3D_Real X_Length,Y_Length,Z_Length,Volume
DIF3D_Real :: HexC1 = 0.5D0, HexC2 = 0.866025403784439d0, HexC3 = 0.9306048591021d0
DIF3D_Real HexPitch


100 FORMAT('[GEODST]...SORRY, BUT I MUST STOP')
101 FORMAT('[GEODST]',107('.'))
102 FORMAT('')
105 FORMAT('[GEODST]...There was a fatal error that occured in (ComputeAverageThickness)',56('.'))

IF (.NOT. User_GEODST%Defined) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...GEODST structure must be defined in order for this to work :-/ ")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

CALL GEODST_Copy(User_GEODST,Copy_GEODST)
CALL GEODST_NRASS1(Copy_GEODST)

ALLOCATE(Volumes(Copy_GEODST%NREG),STAT=IOS)
IF (IOS .NE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...Scratch memory failure 1")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF
Volumes = 0.0D0
Thickness = 0.0d0

DO K = 1,Copy_GEODST%NINTK
   DO J = 1,Copy_GEODST%NINTJ
      DO I = 1,Copy_GEODST%NINTI
         III = Copy_GEODST%MR(I,J,K) ! The region number assigned to this mesh
         IF (III .NE. 0) THEN
            IF ((Copy_GEODST%IGOM .EQ. 1) .OR. (Copy_GEODST%IGOM .EQ. 6) .OR. (Copy_GEODST%IGOM .EQ. 14)) THEN  ! Generic cartesian x-y-z mesh option
               X_Length = Copy_GEODST%XMESH(I+1) - Copy_GEODST%XMESH(I)
               Y_Length = Copy_GEODST%YMESH(J+1) - Copy_GEODST%YMESH(J)
               Z_Length = Copy_GEODST%ZMESH(K+1) - Copy_GEODST%ZMESH(K)
               IF (Copy_GEODST%IGOM .EQ. 1) Y_Length = 1.0d0
               IF (Copy_GEODST%IGOM .LE. 6) Z_Length = 1.0d0
               Volume = X_Length*Y_Length*Z_Length
               Thickness(1,III) = Thickness(1,III) + X_Length*Volume
               Thickness(2,III) = Thickness(2,III) + Y_Length*Volume
               Thickness(3,III) = Thickness(3,III) + Z_Length*Volume
               Volumes(III) = Volumes(III) + Volume
            ELSE IF ((Copy_GEODST%IGOM .EQ. 10) .OR. (Copy_GEODST%IGOM .EQ. 18)) THEN ! HexZ 
               HexPitch = Copy_GEODST%XMESH(2)
               X_Length = HexPitch * HexC3
               Y_Length = X_Length
               Z_Length = Copy_GEODST%ZMESH(K+1) - Copy_GEODST%ZMESH(K)
               IF (Copy_GEODST%IGOM .LE. 10) Z_Length = 1.0d0
               Volume = X_Length*Y_Length*Z_Length
               Thickness(1,III) = Thickness(1,III) + X_Length * Volume
               Thickness(2,III) = Thickness(2,III) + Y_Length * Volume
               Thickness(3,III) = Thickness(3,III) + Z_Length * Volume
               Volumes(III) = Volumes(III) + Volume
            ELSE IF (Copy_GEODST%IGOM .EQ. 17) THEN ! TriZ
               IOS = (I+1)/2   !  1,2,3,4,5,6,.... -> 1,1,2,2,3,3,4,4,...   * 2 = 2,2,4,4,6,6,8,8,...
               ! Use these as pitch
               X_Length = 0.288675134594813d0 * HexPitch
               Y_Length = 0.5d0 * HexPitch
               Z_Length = Copy_GEODST%ZMESH(K+1) - Copy_GEODST%ZMESH(K)
               Volume = X_Length*Y_Length*Z_Length
               Thickness(1,III) = Thickness(1,III) + X_Length * Volume
               Thickness(2,III) = Thickness(2,III) + Y_Length * Volume
               Thickness(3,III) = Thickness(3,III) + Z_Length * Volume
               Volumes(III) = Volumes(III) + Volume
            ELSE IF (Copy_GEODST%IGOM .EQ. 7) THEN ! RZ
               X_Length = Copy_GEODST%XMESH(I+1) - Copy_GEODST%XMESH(I)
               Y_Length = Copy_GEODST%YMESH(J+1) - Copy_GEODST%YMESH(J)
               Z_Length = 1.0d0 !Copy_GEODST%ZMESH(K+1) - Copy_GEODST%ZMESH(K)
               Volume = (Copy_GEODST%XMESH(I+1)*Copy_GEODST%XMESH(I+1) &
                      -  Copy_GEODST%XMESH(I)*Copy_GEODST%XMESH(I)     )*DIF3D_Definition_PI*Y_Length*Z_Length
               !WRITE(MODULE_OUT,'(" X Y Z =",3(1PE13.6,1X)," Volume ",1PE13.6)') &
               !   X_Length,Y_Length,Z_Length,Volume
               Thickness(1,III) = Thickness(1,III) + X_Length * Volume
               Thickness(2,III) = Thickness(2,III) + Y_Length * Volume
               Thickness(3,III) = Thickness(3,III) + Z_Length * Volume
               Volumes(III) = Volumes(III) + Volume
            ELSE
               WRITE(MODULE_OUT,105)
               WRITE(MODULE_OUT,'("[GEODST]...Unknown geometry type, Call MAS :-( ")')
               WRITE(MODULE_OUT,100)
               CALL Abort
            END IF
         END IF
      END DO
   END DO
END DO

!WRITE(MODULE_OUT,*)'GEODST: AverageThickness'
!DO I = 1,Copy_GEODST%NREG
!   WRITE(MODULE_OUT,'("Region ",I5," thickness ",3(1PE13.6,1X)," Volume ",1PE13.6)') &
!      I,Thickness(1,I),Thickness(2,I),Thickness(3,I),Volume
!END DO

! To finish the centroid, we must divide by the volume
DO I = 1,Copy_GEODST%NREG
   Volume = Volumes(I)
   IF (Volume .NE. 0.0d0) Volume = 1.0d0/Volume
   Thickness(1,I) = Thickness(1,I) * Volume
   Thickness(2,I) = Thickness(2,I) * Volume
   Thickness(3,I) = Thickness(3,I) * Volume
END DO

DEALLOCATE(Volumes,STAT=IOS)
IF (IOS .NE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...Scratch memory failure 2")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

CALL GEODST_Void(Copy_GEODST)

END SUBROUTINE GEODST_ComputeAverageThickness

!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_MeshVolumes(User_GEODST,MeshVolumes)
!  This routine computes and returns the mesh by mesh volume in the GEODST mesh
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GEODST_MeshVolumes(User_GEODST,MeshVolumes)
IMPLICIT NONE
! Passed in
TYPE (GEODST_DATA) User_GEODST                  ! A user GEODST data structure
DIF3D_Real        :: MeshVolumes(User_GEODST%NINTI,User_GEODST%NINTJ,User_GEODST%NINTK)
! Local junk
DIF3D_Int I,J,K,III,IOS
DIF3D_Real X_Length,Y_Length,Z_Length,Volume
DIF3D_Real :: HexC1 = 0.5D0, HexC2 = 0.866025403784439d0, HexC3 = 0.9306048591021d0
DIF3D_Real HexPitch

100 FORMAT('[GEODST]...SORRY, BUT I MUST STOP')
101 FORMAT('[GEODST]',107('.'))
102 FORMAT('')
105 FORMAT('[GEODST]...There was a fatal error that occured in (MeshVolumes)',56('.'))

IF (.NOT. User_GEODST%Defined) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...GEODST structure must be defined in order for this to work :-/ ")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (User_GEODST%NRASS .NE. 1) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...NRASS must be equal to 1 in this subroutine :-/ ")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

MeshVolumes = 0.0d0

DO K = 1,User_GEODST%NINTK
   DO J = 1,User_GEODST%NINTJ
      DO I = 1,User_GEODST%NINTI
         III = User_GEODST%MR(I,J,K) ! The region number assigned to this mesh
         IF (III .NE. 0) THEN
            IF ((User_GEODST%IGOM .EQ. 1) .OR. (User_GEODST%IGOM .EQ. 6) .OR. (User_GEODST%IGOM .EQ. 14)) THEN  ! Generic cartesian x-y-z mesh option
               X_Length = User_GEODST%XMESH(I+1) - User_GEODST%XMESH(I)
               Y_Length = User_GEODST%YMESH(J+1) - User_GEODST%YMESH(J)
               Z_Length = User_GEODST%ZMESH(K+1) - User_GEODST%ZMESH(K)
               IF (User_GEODST%IGOM .EQ. 1) Y_Length = 1.0d0
               IF (User_GEODST%IGOM .LE. 6) Z_Length = 1.0d0
               MeshVolumes(I,J,K) = X_Length*Y_Length*Z_Length
            ELSE IF ((User_GEODST%IGOM .EQ. 10) .OR. (User_GEODST%IGOM .EQ. 18)) THEN ! HexZ
               HexPitch = User_GEODST%XMESH(2)
               X_Length = HexPitch * HexC3
               Y_Length = X_Length
               Z_Length = User_GEODST%ZMESH(K+1) - User_GEODST%ZMESH(K)
               IF (User_GEODST%IGOM .EQ. 10) Z_Length = 1.0d0
               MeshVolumes(I,J,K) = X_Length*Y_Length*Z_Length
            ELSE IF (User_GEODST%IGOM .EQ. 17) THEN ! TriZ
               HexPitch = User_GEODST%XMESH(2)
               X_Length = 0.288675134594813d0 * HexPitch
               Y_Length = 0.5d0 * HexPitch
               Z_Length = User_GEODST%ZMESH(K+1) - User_GEODST%ZMESH(K)
               MeshVolumes(I,J,K) = X_Length*Y_Length*Z_Length
            ELSE IF (User_GEODST%IGOM .EQ. 7) THEN ! RZ
               X_Length = User_GEODST%XMESH(I+1) - User_GEODST%XMESH(I)
               Y_Length = User_GEODST%YMESH(J+1) - User_GEODST%YMESH(J)
               Z_Length = 1.0d0 ! User_GEODST%ZMESH(K+1) - User_GEODST%ZMESH(K)
               MeshVolumes(I,J,K) = (User_GEODST%XMESH(I+1)*User_GEODST%XMESH(I+1) &
                                  -  User_GEODST%XMESH(I)*User_GEODST%XMESH(I)     )*DIF3D_Definition_PI*Y_Length*Z_Length
            ELSE IF (User_GEODST%IGOM .EQ. 15) THEN ! ThetaRZ = Area*Z = pi*R^2 * Z = Delta_Theta/2 * (R2^2 - R1^2)*Z
               X_Length = (User_GEODST%XMESH(I+1) - User_GEODST%XMESH(I))*0.5D0 ! Theta (0,2PI)
               Y_Length = User_GEODST%YMESH(J+1) - User_GEODST%YMESH(J) ! R
               Z_Length = User_GEODST%ZMESH(K+1) - User_GEODST%ZMESH(K) ! Z
               MeshVolumes(I,J,K) = (User_GEODST%YMESH(J+1)*User_GEODST%YMESH(J+1) &
                                  -  User_GEODST%YMESH(J)*User_GEODST%YMESH(J)     )*X_Length*Z_Length
            ELSE
               WRITE(MODULE_OUT,105)
               WRITE(MODULE_OUT,'("[GEODST]...Unknown geometry type, Call MAS :-( ")')
               WRITE(MODULE_OUT,100)
               CALL Abort
            END IF
         END IF
      END DO
   END DO
END DO

END SUBROUTINE GEODST_MeshVolumes

!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_RegionsOnBoundary(User_GEODST,Coordinates)
!  This routine computes and returns a vector that defines which regions touch which boundaries and roughly at which position
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GEODST_RegionsOnBoundary(User_GEODST,Coordinates)
IMPLICIT NONE
! Passed in
TYPE (GEODST_DATA) User_GEODST                  ! A user GEODST data structure
! Local arrays
DIF3D_Real Coordinates(3,6,User_GEODST%NREG)    ! 1:3 are the center coordinates of the surface
                                                ! 1:3 are for +/-X, +/-Y, +/-Z
DIF3D_Real, POINTER :: SurfacesPerRegion(:,:)
TYPE (GEODST_DATA) Copy_GEODST                  ! A user GEODST data structure
! Local junk
DIF3D_Int I,J,K,iRegion,III,IOS
DIF3D_Real :: HexC1 = 0.5D0, HexC2 = 0.866025403784439d0, HexC3 = 0.9306048591021d0
DIF3D_Real HexPitch
!  RZ                         Cartesian                   Hex                    Tri-Z
!  -X does not exist         +/- X is obvious                +Y                  same idea as Hex
!  +X last radial ring       +/- Y is obvious                __
!  -Y lower boundary         +/- Z is obvious            -X /  \ +X
!  +Y upper boundary                                        \__/
!                                                            -Y
!                                                       +/- are obvious       

100 FORMAT('[GEODST]...SORRY, BUT I MUST STOP')
101 FORMAT('[GEODST]',107('.'))
102 FORMAT('')
105 FORMAT('[GEODST]...There was a fatal error that occured in (RegionsOnBoundary)',56('.'))

IF (.NOT. User_GEODST%Defined) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...GEODST structure must be defined in order for this to work :-/ ")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

CALL GEODST_Copy(User_GEODST,Copy_GEODST)
CALL GEODST_NRASS1(Copy_GEODST)

ALLOCATE(SurfacesPerRegion(6,Copy_GEODST%NREG),STAT=IOS)
IF (IOS .NE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...Scratch memory failure 1")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

Coordinates = 0.0d0
SurfacesPerRegion = 0.0d0
IF ((Copy_GEODST%IGOM .EQ. 1) .OR. & ! x
    (Copy_GEODST%IGOM .EQ. 6) .OR. & ! x-y
    (Copy_GEODST%IGOM .EQ. 7) .OR. & ! RZ
    (Copy_GEODST%IGOM .EQ. 14)     ) THEN ! x-y-z
   ! -X boundary
   DO J = 1,Copy_GEODST%NINTJ
      IF ((Copy_GEODST%IMB1 .EQ. 0) .OR. (Copy_GEODST%IMB1 .EQ. 2)) THEN
      III = 0
      DO I = 1,Copy_GEODST%NINTI
         IF (Copy_GEODST%MR(I,J,1) .NE. 0) THEN
            III = I
            EXIT
         END IF
      END DO
      IF (III .NE. 0) THEN
         DO K = 1,Copy_GEODST%NINTK
            iRegion = Copy_GEODST%MR(I,J,K)
            SurfacesPerRegion(1,iRegion) = SurfacesPerRegion(1,iRegion) + 1.0d0
            Coordinates(1,1,iRegion) = Coordinates(1,1,iRegion) + Copy_GEODST%XMESH(I)
            Coordinates(2,1,iRegion) = Coordinates(2,1,iRegion) + 0.5D0*(Copy_GEODST%YMESH(J+1)+Copy_GEODST%YMESH(J))
            Coordinates(3,1,iRegion) = Coordinates(3,1,iRegion) + 0.5D0*(Copy_GEODST%ZMESH(K+1)+Copy_GEODST%ZMESH(K))
         END DO
      END IF
      END IF
   !END DO
   ! +X boundary
   !DO J = 1,Copy_GEODST%NINTJ
      IF ((Copy_GEODST%IMB2 .EQ. 0) .OR. (Copy_GEODST%IMB2 .EQ. 2)) THEN
      III = 0
      DO I = Copy_GEODST%NINTI,1,-1
         IF (Copy_GEODST%MR(I,J,1) .NE. 0) THEN
            III = I
            EXIT
         END IF
      END DO
      IF (III .NE. 0) THEN
         DO K = 1,Copy_GEODST%NINTK
            iRegion = Copy_GEODST%MR(I,J,K)
            SurfacesPerRegion(2,iRegion) = SurfacesPerRegion(2,iRegion) + 1.0d0
            Coordinates(1,2,iRegion) = Coordinates(1,2,iRegion) + Copy_GEODST%XMESH(I+1)
            Coordinates(2,2,iRegion) = Coordinates(2,2,iRegion) + 0.5D0*(Copy_GEODST%YMESH(J+1)+Copy_GEODST%YMESH(J))
            Coordinates(3,2,iRegion) = Coordinates(3,2,iRegion) + 0.5D0*(Copy_GEODST%ZMESH(K+1)+Copy_GEODST%ZMESH(K))
         END DO
      END IF
      END IF
   END DO
   ! -Y boundary
   DO I = 1,Copy_GEODST%NINTI
      IF ((Copy_GEODST%JMB1 .EQ. 0) .OR. (Copy_GEODST%JMB1 .EQ. 2)) THEN
      III = 0
      DO J = 1,Copy_GEODST%NINTJ
         IF (Copy_GEODST%MR(I,J,1) .NE. 0) THEN
            III = J
            EXIT
         END IF
      END DO
      IF (III .NE. 0) THEN
         DO K = 1,Copy_GEODST%NINTK
            iRegion = Copy_GEODST%MR(I,J,K)
            SurfacesPerRegion(3,iRegion) = SurfacesPerRegion(3,iRegion) + 1.0d0
            Coordinates(1,3,iRegion) = Coordinates(1,3,iRegion) + 0.5D0*(Copy_GEODST%XMESH(I+1)+Copy_GEODST%XMESH(I))
            Coordinates(2,3,iRegion) = Coordinates(2,3,iRegion) + Copy_GEODST%YMESH(J)
            Coordinates(3,3,iRegion) = Coordinates(3,3,iRegion) + 0.5D0*(Copy_GEODST%ZMESH(K+1)+Copy_GEODST%ZMESH(K))
         END DO
      END IF
      END IF
   !END DO
   ! +Y boundary
   !DO I = 1,Copy_GEODST%NINTI
      IF ((Copy_GEODST%JMB2 .EQ. 0) .OR. (Copy_GEODST%JMB2 .EQ. 2)) THEN
      III = 0
      DO J = Copy_GEODST%NINTJ,1,-1
         IF (Copy_GEODST%MR(I,J,1) .NE. 0) THEN
            III = J
            EXIT
         END IF
      END DO
      IF (III .NE. 0) THEN
         DO K = 1,Copy_GEODST%NINTK
            iRegion = Copy_GEODST%MR(I,J,K)
            SurfacesPerRegion(4,iRegion) = SurfacesPerRegion(4,iRegion) + 1.0d0
            Coordinates(1,4,iRegion) = Coordinates(1,4,iRegion) + 0.5D0*(Copy_GEODST%XMESH(I+1)+Copy_GEODST%XMESH(I))
            Coordinates(2,4,iRegion) = Coordinates(2,4,iRegion) + Copy_GEODST%YMESH(J+1)
            Coordinates(3,4,iRegion) = Coordinates(3,4,iRegion) + 0.5D0*(Copy_GEODST%ZMESH(K+1)+Copy_GEODST%ZMESH(K))
         END DO
      END IF
      END IF
   END DO
   IF (Copy_GEODST%IGOM .EQ. 14) THEN 
      ! +/- Z boundary
      DO I = 1,Copy_GEODST%NINTI
         DO J = 1,Copy_GEODST%NINTJ
            K = 1
            IF ((Copy_GEODST%MR(I,J,K) .NE. 0) .AND. ((Copy_GEODST%KMB1 .EQ. 0) .OR. (Copy_GEODST%KMB1 .EQ. 2)) ) THEN
               iRegion = Copy_GEODST%MR(I,J,K)
               SurfacesPerRegion(5,iRegion) = SurfacesPerRegion(5,iRegion) + 1.0d0
               Coordinates(1,5,iRegion) = Coordinates(1,5,iRegion) + 0.5D0*(Copy_GEODST%XMESH(I+1)+Copy_GEODST%XMESH(I))
               Coordinates(2,5,iRegion) = Coordinates(2,5,iRegion) + 0.5D0*(Copy_GEODST%YMESH(J+1)+Copy_GEODST%YMESH(J))
               Coordinates(3,5,iRegion) = Coordinates(3,5,iRegion) + Copy_GEODST%ZMESH(K)
            END IF
            K = Copy_GEODST%NINTK
            IF ((Copy_GEODST%MR(I,J,K) .NE. 0) .AND. ((Copy_GEODST%KMB2 .EQ. 0) .OR. (Copy_GEODST%KMB2 .EQ. 2)) ) THEN
               iRegion = Copy_GEODST%MR(I,J,K)
               SurfacesPerRegion(6,iRegion) = SurfacesPerRegion(6,iRegion) + 1.0d0
               Coordinates(1,6,iRegion) = Coordinates(1,6,iRegion) + 0.5D0*(Copy_GEODST%XMESH(I+1)+Copy_GEODST%XMESH(I))
               Coordinates(2,6,iRegion) = Coordinates(2,6,iRegion) + 0.5D0*(Copy_GEODST%YMESH(J+1)+Copy_GEODST%YMESH(J))
               Coordinates(3,6,iRegion) = Coordinates(3,6,iRegion) + Copy_GEODST%ZMESH(K+1)
            END IF
         END DO
      END DO
   END IF
ELSE IF ((Copy_GEODST%IGOM .EQ. 10) .OR. (Copy_GEODST%IGOM .EQ. 18)) THEN ! HexZ
   HexPitch = Copy_GEODST%XMESH(2)
   ! -X boundary
   DO J = 1,Copy_GEODST%NINTJ
      III = 0
      DO I = 1,Copy_GEODST%NINTI
         IF (Copy_GEODST%MR(I,J,1) .NE. 0) THEN
            III = I
            EXIT
         END IF
      END DO
      IF (III .NE. 0) THEN
         DO K = 1,Copy_GEODST%NINTK
            iRegion = Copy_GEODST%MR(I,J,K)
            SurfacesPerRegion(1,iRegion) = SurfacesPerRegion(1,iRegion) + 1.0d0
            Coordinates(1,1,iRegion) = Coordinates(1,1,iRegion) + ((I-1)*HexPitch - (J-1)*HexC1*HexPitch) - 0.5d0*HexPitch
            Coordinates(2,1,iRegion) = Coordinates(2,1,iRegion) + (                 (J-1)*HexC2*HexPitch)
            Coordinates(3,1,iRegion) = Coordinates(3,1,iRegion) + 0.5d0 * (Copy_GEODST%ZMESH(K+1) - Copy_GEODST%ZMESH(K))
         END DO
      END IF
   !END DO
   ! +X boundary
   !DO J = 1,Copy_GEODST%NINTJ
      III = 0
      DO I = Copy_GEODST%NINTI,1,-1
         IF (Copy_GEODST%MR(I,J,1) .NE. 0) THEN
            III = I
            EXIT
         END IF
      END DO
      IF (III .NE. 0) THEN
         DO K = 1,Copy_GEODST%NINTK
            iRegion = Copy_GEODST%MR(I,J,K)
            SurfacesPerRegion(2,iRegion) = SurfacesPerRegion(2,iRegion) + 1.0d0
            Coordinates(1,2,iRegion) = Coordinates(1,2,iRegion) + ((I-1)*HexPitch - (J-1)*HexC1*HexPitch) + 0.5d0*HexPitch
            Coordinates(2,2,iRegion) = Coordinates(2,2,iRegion) + (                 (J-1)*HexC2*HexPitch)
            Coordinates(3,2,iRegion) = Coordinates(3,2,iRegion) + 0.5d0 * (Copy_GEODST%ZMESH(K+1) - Copy_GEODST%ZMESH(K))
         END DO
      END IF
   END DO
   ! -Y boundary
   DO I = 1,Copy_GEODST%NINTI
      III = 0
      DO J = 1,Copy_GEODST%NINTJ
         IF (Copy_GEODST%MR(I,J,1) .NE. 0) THEN
            III = J
            EXIT
         END IF
      END DO
      IF (III .NE. 0) THEN
         DO K = 1,Copy_GEODST%NINTK
            iRegion = Copy_GEODST%MR(I,J,K)
            SurfacesPerRegion(3,iRegion) = SurfacesPerRegion(3,iRegion) + 1.0d0
            Coordinates(1,3,iRegion) = Coordinates(1,3,iRegion) + ((I-1)*HexPitch - (J-1)*HexC1*HexPitch) 
            Coordinates(2,3,iRegion) = Coordinates(2,3,iRegion) + (                 (J-1)*HexC2*HexPitch) - 0.5d0*HexC2*HexPitch
            Coordinates(3,3,iRegion) = Coordinates(3,3,iRegion) + 0.5d0 * (Copy_GEODST%ZMESH(K+1) - Copy_GEODST%ZMESH(K))
         END DO
      END IF
   !END DO
   ! +Y boundary
   !DO I = 1,Copy_GEODST%NINTI
      III = 0
      DO J = Copy_GEODST%NINTJ,1,-1
         IF (Copy_GEODST%MR(I,J,1) .NE. 0) THEN
            III = J
            EXIT
         END IF
      END DO
      IF (III .NE. 0) THEN
         DO K = 1,Copy_GEODST%NINTK
            iRegion = Copy_GEODST%MR(I,J,K)
            SurfacesPerRegion(4,iRegion) = SurfacesPerRegion(4,iRegion) + 1.0d0
            Coordinates(1,4,iRegion) = Coordinates(1,4,iRegion) + ((I-1)*HexPitch - (J-1)*HexC1*HexPitch)
            Coordinates(2,4,iRegion) = Coordinates(2,4,iRegion) + (                 (J-1)*HexC2*HexPitch) + 0.5d0*HexC2*HexPitch
            Coordinates(3,4,iRegion) = Coordinates(3,4,iRegion) + 0.5d0 * (Copy_GEODST%ZMESH(K+1) - Copy_GEODST%ZMESH(K))
         END DO
      END IF
   END DO
   ! +/- Z boundary
   DO I = 1,Copy_GEODST%NINTI
      DO J = 1,Copy_GEODST%NINTJ
         K = 1
         IF ((Copy_GEODST%MR(I,J,K) .NE. 0) .AND. ((Copy_GEODST%KMB1 .EQ. 0) .OR. (Copy_GEODST%KMB1 .EQ. 2)) ) THEN
            iRegion = Copy_GEODST%MR(I,J,K)
            SurfacesPerRegion(5,iRegion) = SurfacesPerRegion(5,iRegion) + 1.0d0
            Coordinates(1,5,iRegion) = Coordinates(1,5,iRegion) + ((I-1)*HexPitch - (J-1)*HexC1*HexPitch)
            Coordinates(2,5,iRegion) = Coordinates(2,5,iRegion) + (                 (J-1)*HexC2*HexPitch)
            Coordinates(3,5,iRegion) = Coordinates(3,5,iRegion) + Copy_GEODST%ZMESH(K)
         END IF
         K = Copy_GEODST%NINTK
         IF ((Copy_GEODST%MR(I,J,K) .NE. 0) .AND. ((Copy_GEODST%KMB2 .EQ. 0) .OR. (Copy_GEODST%KMB2 .EQ. 2)) ) THEN
            iRegion = Copy_GEODST%MR(I,J,K)
            SurfacesPerRegion(6,iRegion) = SurfacesPerRegion(6,iRegion) + 1.0d0
            Coordinates(1,6,iRegion) = Coordinates(1,6,iRegion) + ((I-1)*HexPitch - (J-1)*HexC1*HexPitch)
            Coordinates(2,6,iRegion) = Coordinates(2,6,iRegion) + (                 (J-1)*HexC2*HexPitch)
            Coordinates(3,6,iRegion) = Coordinates(3,6,iRegion) + Copy_GEODST%ZMESH(K+1)
         END IF
      END DO
   END DO
!ELSE IF (Copy_GEODST%IGOM .EQ. 17) THEN ! TriZ
ELSE
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...Unknown geometry type, Call MAS :-( ")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

!DO I = 1,Copy_GEODST%NREG
!   DO J = 1,6
!      WRITE(MODULE_OUT,'("[GEODST]...Region ",I5," boundary ",I5," touched ",I5, &
!                         " surfaces which have summed coordinates ",3(1PE13.6,1X))') &
!      I,J,SurfacesPerRegion(J,I),(Coordinates(K,J,I),K=1,3)
!   END DO
!END DO

! Compute the average surface of each mesh point
DO I = 1,Copy_GEODST%NREG
   DO J = 1,6
      IF (SurfacesPerRegion(J,I) .NE. 0.0d0) THEN
         Coordinates(1,J,I) = Coordinates(1,J,I)/SurfacesPerRegion(J,I)
         Coordinates(2,J,I) = Coordinates(2,J,I)/SurfacesPerRegion(J,I)
         Coordinates(3,J,I) = Coordinates(3,J,I)/SurfacesPerRegion(J,I)
      ELSE
         Coordinates(1,J,I) = NE_FreeForm_Error_Real
         Coordinates(2,J,I) = NE_FreeForm_Error_Real
         Coordinates(3,J,I) = NE_FreeForm_Error_Real
      END IF
   END DO
END DO

DEALLOCATE(SurfacesPerRegion,STAT=IOS)
IF (IOS .NE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...Scratch memory failure 2")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

CALL GEODST_Void(Copy_GEODST)

END SUBROUTINE GEODST_RegionsOnBoundary

!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_NODAL_factors(User_GEODST,MeshFluxFactors)
!  This routine computes and returns the mesh by mesh NODAL flux factors
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GEODST_NODAL_factors(User_GEODST,MeshFluxFactors)
IMPLICIT NONE
! Passed in
TYPE (GEODST_DATA) User_GEODST                  ! A user GEODST data structure
DIF3D_Real        :: MeshFluxFactors(User_GEODST%NINTI,User_GEODST%NINTJ)
! Local junk
DIF3D_Int I,J,HexOption
DIF3D_Real Delta

100 FORMAT('[GEODST]...SORRY, BUT I MUST STOP')
101 FORMAT('[GEODST]',107('.'))
102 FORMAT('')
105 FORMAT('[GEODST]...There was a fatal error that occured in (MeshVolumes)',56('.'))

IF (.NOT. User_GEODST%Defined) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[GEODST]...GEODST structure must be defined in order for this to work :-/ ")')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

MeshFluxFactors = 1.0d0 ! All other geometry cases

! I need to detect whether this is 
! HexOption = 0 is a full core
! HexOption = 1 is 120 symmetry/periodic which will have Local_GEODST%IMB1=4 and Local_GEODST%NTRIAG = 0
! HexOption = 2 is  60 symmetry/periodic which will have Local_GEODST%IMB1=4 and Local_GEODST%NTRIAG = 1
HexOption = 0
IF ((User_GEODST%IGOM .EQ. 10) .OR. (User_GEODST%IGOM .EQ. 18)) THEN
   IF (User_GEODST%IMB1 .EQ. 4) THEN
      HexOption = 1
      IF (User_GEODST%NTRIAG .EQ. 1) HexOption = 2
   END IF
END IF
IF ((HexOption .NE. 0) .AND. (User_GEODST%NINTJ .GT. 1)) THEN
   ! We translate the I=1 J=2:NINTJ part as those GEODST cells will not be used
   DO I = 1,User_GEODST%NINTI
      DO J = 1,User_GEODST%NINTJ
         MeshFluxFactors(I,J) = 1.0d0
      END DO
   END DO
   Delta = 3.0d0 ! 120 sym/per
   IF (HexOption .EQ. 2) Delta = 6.0d0 ! 60 sym/per
   MeshFluxFactors(1,1) = 1.0d0/Delta
   ! Correct the nodal flux map and eliminate their contribution
   DO J = 2,User_GEODST%NINTJ
      MeshFluxFactors(1,J) = 0.0d0
   END DO
END IF

END SUBROUTINE GEODST_NODAL_factors

!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_VOID(User_GEODST)
!  Provides a path to void a GEODST data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GEODST_VOID(User_GEODST)
IMPLICIT NONE
TYPE (GEODST_DATA) User_GEODST           ! A user data variable to be defined by reading in a GEODST file
! local
DIF3D_Int IOS
   
100 FORMAT('[GEODST]...SORRY, BUT I MUST STOP')
101 FORMAT('[GEODST]',107('.'))
102 FORMAT('')
105 FORMAT('[GEODST]...There was a fatal error that occured in (VOID_GEODSTTYPE)',31('.'),A16)
   
IF (User_GEODST%DEFINED) THEN
   DEALLOCATE(User_GEODST%XMESH,User_GEODST%IFINTS,User_GEODST%YMESH,User_GEODST%JFINTS,                        &
              User_GEODST%ZMESH,User_GEODST%KFINTS,User_GEODST%VOLR,User_GEODST%BSQ,User_GEODST%BNDC,User_GEODST%BNCI,  &
              User_GEODST%NZHBB,User_GEODST%NZC,User_GEODST%NZNR,STAT = IOS                                 )
   IF (IOS .NE. 0) THEN
      WRITE(MODULE_OUT,105)
      WRITE(MODULE_OUT,'("[GEODST]...FAILED TO DEALLOCATE MEMORY FOR GEODST",66("."))')
      WRITE(MODULE_OUT,100)
      CALL Abort
   END IF 
   DEALLOCATE(User_GEODST%MR,STAT = IOS)
   ! Ignore the error
END IF 
   
User_GEODST%DEFINED  = .FALSE.
END SUBROUTINE GEODST_VOID

!---------------------------------------------------------------------------------------------------------------------------------
!  GEODST_PRINT(User_GEODST)
!  Prints the entire GEODST data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GEODST_PRINT(User_GEODST)
IMPLICIT NONE
! Passed in
TYPE (GEODST_DATA) User_GEODST                  ! A user GEODST data structure
! Local junk
DIF3D_Int IOS,I,J,K,N

100 FORMAT('[GEODST]...SORRY, BUT I MUST STOP')
101 FORMAT('[GEODST]',107('.'))
102 FORMAT('')
105 FORMAT('[GEODST]...There was a fatal error that occured in (PRINT)',57('.'))

IF (.NOT. User_GEODST%DEFINED) THEN
   WRITE(MODULE_OUT,101)
   WRITE(MODULE_OUT,'("[GEODST]...GEODST DATA STRCTURE NOT DEFINED IN PRINT...NON-FATAL",35("."),A16)') User_GEODST%FILENAME
   WRITE(MODULE_OUT,101)
ELSE
   WRITE(MODULE_OUT,'("[GEODST]...HNAME...........................................",40("."),A16)') User_GEODST%HNAME
   WRITE(MODULE_OUT,'("[GEODST]...HUSE(1).........................................",40("."),A16)') User_GEODST%HUSE(1)
   WRITE(MODULE_OUT,'("[GEODST]...HUSE(2).........................................",40("."),A16)') User_GEODST%HUSE(2)
   WRITE(MODULE_OUT,'("[GEODST]...IVERS...........................................",40("."),I16)') User_GEODST%IVERS
   WRITE(MODULE_OUT,'("[GEODST]...GEOMETRY TYPE...................................",40("."),I16)') User_GEODST%IGOM
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF ZONES.................................",40("."),I16)') User_GEODST%NZONE
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF REGIONS...............................",40("."),I16)') User_GEODST%NREG
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF ZONE CLASSIFICATIONS..................",40("."),I16)') User_GEODST%NZCL
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF FIRST  DIMENSION COARSE MESH INTERVALS",40("."),I16)') User_GEODST%NCINTI
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF SECOND DIMENSION COARSE MESH INTERVALS",40("."),I16)') User_GEODST%NCINTJ
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF THIRD  DIMENSION COARSE MESH INTERVALS",40("."),I16)') User_GEODST%NCINTK
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF FIRST  DIMENSION FINE   MESH INTERVALS",40("."),I16)') User_GEODST%NINTI
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF SECOND DIMENSION FINE   MESH INTERVALS",40("."),I16)') User_GEODST%NINTJ
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF THIRD  DIMENSION FINE   MESH INTERVALS",40("."),I16)') User_GEODST%NINTK
   WRITE(MODULE_OUT,'("[GEODST]...FIRST  BOUNDARY CONDITION ON FIRST  DIMENSION...",40("."),I16)') User_GEODST%IMB1
   WRITE(MODULE_OUT,'("[GEODST]...SECOND BOUNDARY CONDITION ON FIRST  DIMENSION...",40("."),I16)') User_GEODST%IMB2
   WRITE(MODULE_OUT,'("[GEODST]...FIRST  BOUNDARY CONDITION ON SECOND DIMENSION...",40("."),I16)') User_GEODST%JMB1
   WRITE(MODULE_OUT,'("[GEODST]...SECOND BOUNDARY CONDITION ON SECOND DIMENSION...",40("."),I16)') User_GEODST%JMB2
   WRITE(MODULE_OUT,'("[GEODST]...FIRST  BOUNDARY CONDITION ON THIRD  DIMENSION...",40("."),I16)') User_GEODST%KMB1
   WRITE(MODULE_OUT,'("[GEODST]...SECOND BOUNDARY CONDITION ON THIRD  DIMENSION...",40("."),I16)') User_GEODST%KMB2
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF BUCKLING SPECIFICATIONS...............",40("."),I16)') User_GEODST%NBS
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF CONSTANTS FOR EXTERNAL BOUNDARIES.....",40("."),I16)') User_GEODST%NBCS
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF CONSTANTS FOR INTERNAL BOUNDARIES.....",40("."),I16)') User_GEODST%NIBCS
   WRITE(MODULE_OUT,'("[GEODST]...NUMBER OF ZONES WHICH ARE BLACK ABSORBERS.......",40("."),I16)') User_GEODST%NZWBB
   WRITE(MODULE_OUT,'("[GEODST]...TRIANGULAR/HEXAGONAL GEOMETRY OPTION............",40("."),I16)') User_GEODST%NTRIAG
   WRITE(MODULE_OUT,'("[GEODST]...REGION ASSIGNMENT...............................",40("."),I16)') User_GEODST%NRASS
   WRITE(MODULE_OUT,'("[GEODST]...ORIENTATION OF FIRST FINE MESH INTERVAL (TRI)...",40("."),I16)') User_GEODST%NTHPT
   WRITE(MODULE_OUT,'("[GEODST]...RESERVED 1......................................",40("."),I16)') User_GEODST%NGOP(1)
   WRITE(MODULE_OUT,'("[GEODST]...RESERVED 2......................................",40("."),I16)') User_GEODST%NGOP(2)
   WRITE(MODULE_OUT,'("[GEODST]...RESERVED 3......................................",40("."),I16)') User_GEODST%NGOP(3)
   WRITE(MODULE_OUT,'("[GEODST]...RESERVED 4......................................",40("."),I16)') User_GEODST%NGOP(4)

   WRITE(MODULE_OUT,'("[GEODST]...X COARSE MESH BOUNDARIES",80("."))')
   DO I = 1,User_GEODST%NCINTI
      WRITE(MODULE_OUT,200) User_GEODST%XMESH(I),User_GEODST%IFINTS(I),User_GEODST%XMESH(I+1)
   END DO
   200 FORMAT('[GEODST]...',F13.6,1X,I3,1X,F13.6)
   !IF (User_GEODST%NCINTJ .GE. 1) THEN
      WRITE(MODULE_OUT,'("[GEODST]...Y COARSE MESH BOUNDARIES",80("."))')
      DO I = 1,User_GEODST%NCINTJ
         WRITE(MODULE_OUT,200) User_GEODST%YMESH(I),User_GEODST%JFINTS(I),User_GEODST%YMESH(I+1)
      END DO
   !END IF
   !IF (User_GEODST%NCINTK .GE. 1) THEN
      WRITE(MODULE_OUT,'("[GEODST]...Z COARSE MESH BOUNDARIES",80("."))')
      DO I = 1,User_GEODST%NCINTK
         WRITE(MODULE_OUT,200) User_GEODST%ZMESH(I),User_GEODST%KFINTS(I),User_GEODST%ZMESH(I+1)
      END DO
   !END IF                     

IF ((User_GEODST%IGOM.GT.0) .OR. (User_GEODST%NBS.GT.0)) THEN
      250 FORMAT('[GEODST]...',5(I4,1X,E13.6,1X))
      255 FORMAT('[GEODST]...',8(I5,1X,I5,1X))
      260 FORMAT('[GEODST]...',11(I5,1X,I2,1X))
      !WRITE(MODULE_OUT,'("[GEODST]...REGION VOLUMES")')
      !WRITE(MODULE_OUT,250) (N,User_GEODST%VOLR(N),N=1,User_GEODST%NREG)
      WRITE(MODULE_OUT,'("[GEODST]...BUCKLING (B^2) VALUES (CM^-2)",75("."))')
      WRITE(MODULE_OUT,250) (N,User_GEODST%BSQ(N),N=1,User_GEODST%NBS)
      WRITE(MODULE_OUT,'("[GEODST]...BOUNDARY CONSTANTS (DEL PHI/PHI =-C/D)",66("."))')
      WRITE(MODULE_OUT,250) (N,User_GEODST%BNDC(N),N=1,User_GEODST%NBCS)
      WRITE(MODULE_OUT,'("[GEODST]...INTERNAL BLACK BOUNDARY CONSTANTS",71("."))')
      WRITE(MODULE_OUT,250) (N,User_GEODST%BNCI(N),N=1,User_GEODST%NIBCS)
      WRITE(MODULE_OUT,'("[GEODST]...ZONE NUMBERS WITH BLACK ABSORBER CONDITIONS",61("."))')
      WRITE(MODULE_OUT,255) (N,User_GEODST%NZHBB(N),N=1,User_GEODST%NZWBB)
      !WRITE(MODULE_OUT,'("[GEODST]...ZONE CLASSIFICATIONS",84("."))')
      !WRITE(MODULE_OUT,260) (N,User_GEODST%NZC(N),N=1,User_GEODST%NZONE)
      !WRITE(MODULE_OUT,'("[GEODST]...ZONE NUMBER ASSIGNED TO EACH REGION",69("."))')
      !WRITE(MODULE_OUT,255) (N,User_GEODST%NZNR(N),N=1,User_GEODST%NREG)
      DO N = 1,User_GEODST%NREG
         I = 0
         IF (User_GEODST%NZNR(N) .NE. 0) I = User_GEODST%NZC(User_GEODST%NZNR(N))
         WRITE(MODULE_OUT,'("[GEODST]...REGION ",I6," IS ASSIGNED ZONE ",I6," OF TYPE ",I5," AND HAS VOLUME ",1PE13.6)') &
            N,User_GEODST%NZNR(N),I,User_GEODST%VOLR(N)
      END DO
END IF

IF ((User_GEODST%IGOM.GT.0) .AND. (User_GEODST%NRASS.EQ.0)) THEN
      !WRITE(MODULE_OUT,'("[GEODST]...REGION NUMBERS Z ASSIGNED TO COARSE MESH INTERVALS I,J,K: (I,J,K)=Z",37("."))')
      300 FORMAT('[GEODST]...',4('(',I3,',',I3,',',I3,')=',I5,3X))
      !WRITE(MODULE_OUT,300) (((I,J,K,User_GEODST%MR(I,J,K),I=1,User_GEODST%NCINTI),J=1,User_GEODST%NCINTJ),K=1,User_GEODST%NCINTK)
      305 FORMAT('[GEODST]...J=',I3,' -> ',3000(1X,I4))
      DO K = 1,User_GEODST%NCINTK
         WRITE(MODULE_OUT,'("[GEODST]...REGION NUMBERS ASSIGNED TO COARSE MESH ON Z PLANE ",I6)') K
         DO J = User_GEODST%NCINTJ,1,-1
            WRITE(MODULE_OUT,305) J,(User_GEODST%MR(I,J,K),I=1,User_GEODST%NCINTI)
         END DO
      END DO
ELSE IF ((User_GEODST%IGOM.GT.0) .AND. (User_GEODST%NRASS.EQ.1)) THEN
      !WRITE(MODULE_OUT,'("[GEODST]...REGION NUMBERS Z ASSIGNED TO FINE MESH INTERVALS I,J,K: (I,J,K)=Z",39("."))')
      !WRITE(MODULE_OUT,300) (((I,J,K,User_GEODST%MR(I,J,K),I=1,User_GEODST%NINTI),J=1,User_GEODST%NINTJ),K=1,User_GEODST%NINTK)
      DO K = 1,User_GEODST%NINTK
         WRITE(MODULE_OUT,'("[GEODST]...REGION NUMBERS ASSIGNED TO COARSE MESH ON Z PLANE ",I6)') K
         DO J = User_GEODST%NINTJ,1,-1
            WRITE(MODULE_OUT,305) J,(User_GEODST%MR(I,J,K),I=1,User_GEODST%NINTI)
         END DO
      END DO
END IF

END IF

END SUBROUTINE GEODST_PRINT

END MODULE GEODST_io

