!---------------------------------------------------------------------------------------------------------------------------------
!  Using the historical format for the RTFLUX file, import/export data from the binary file 
!---------------------------------------------------------------------------------------------------------------------------------
!  RTFLUX_ASSIGNPRINTINFO   !  Sets the output unit and debugprint level settings for the module
!  RTFLUX_DEFINE            !  Allocates and initializes the RTFLUX data structure
!  RTFLUX_VOID              !  Provides a path to deallocate the data
!  RTFLUX_PRINT             !  Prints the RTFLUX data structure
!---------------------------
!  RTFLUX_ASSIGNPRINTINFO(OUTPUT_UNIT)
!  RTFLUX_DEFINE(USER_RTFLUX,NINTI,NINTJ,NINTK,NGROUP)
!  RTFLUX_VOID(USER_RTFLUX)
!  RTFLUX_PRINT(USER_RTFLUX)
!---------------------------------------------------------------------------------------------------------------------------------
MODULE RTFLUX_io
#include "DIF3D_Types.h"
IMPLICIT NONE

TYPE RTFLUX_DATA
   DIF3D_Log       :: DEFINED=.FALSE.
   DIF3D_FileNames       :: FILENAME = 'UNDEFINE'  ! Originating file name
   ! CARD TYPE 0
   DIF3D_Char        :: HNAME    = 'UNDEFINE'
   DIF3D_Char        :: HUSE(2)  = 'UNDEFINE'
   RTFLUX_Int :: IVERS    = 0
   ! CARD TYPE 1
   RTFLUX_Int :: NDIM          !    NUMBER OF DIMENSIONS
   RTFLUX_Int :: NGROUP        !    NUMBER OF ENERGY GROUPS
   RTFLUX_Int :: NINTI         !    NUMBER OF FIRST DIMENSION FINE MESH INTERVALS
   RTFLUX_Int :: NINTJ         !    NUMBER OF SECOND DIMENSION FINE MESH INTERVALS
   RTFLUX_Int :: NINTK         !    NUMBER OF THIRD DIMENSION FINE MESH INTERVALS. NINTK.EQ.1 IF NDIM.LE.2
   RTFLUX_Int :: ITER          !    OUTER ITERATION NUMBER AT WHICH FLUX WAS WRITTEN
   DIF3D_Real_32bit :: EFFK    !    EFFECTIVE MULTIPLICATION FACTOR
   DIF3D_Real_32bit :: POWER   !    POWER IN WATTS TO WHICH FLUX IS NORMALIZED
   RTFLUX_Int :: NBLOK         !    DATA BLOCKING FACTOR
                                       !     IF NDIM.EQ.1 THE GROUP VARIABLE IS BLOCKED INTO NBLOK BLOCKS (SEE 2D RECORD BELOW)
                                       !     IF NDIM.GE.2 THE 2ND DIMENSION VARIABLE IS BLOCKED INTO NBLOK BLOCKS (SEE 3D RECORD)
   ! CARD TYPE 2 & 3
   RTFLUX_Real, DIMENSION(:,:,:,:), POINTER :: FREG  ! (I,J,K,G) REGULAR TOTAL FLUX
END TYPE

TYPE (RTFLUX_DATA), SAVE :: MODULE_RTFLUX ! A optional use "common block" data type for the module
DIF3D_Int            :: MODULE_OUT = 6               ! Module output unit

PRIVATE MODULE_OUT

CONTAINS

!---------------------------------------------------------------------------------------------------------------------------------
!  RTFLUX_ASSIGNPRINTINFO(OUTPUT_UNIT)
!  Sets the output unit and debugprint level settings for the module
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE RTFLUX_ASSIGNPRINTINFO(MODULE_OUTPUT_UNIT)
DIF3D_Int MODULE_OUTPUT_UNIT    ! Reset the object output unit

MODULE_OUT = MODULE_OUTPUT_UNIT

END SUBROUTINE RTFLUX_ASSIGNPRINTINFO

!---------------------------------------------------------------------------------------------------------------------------------
!  RTFLUX_PRINT(USER_RTFLUX)
!  Prints the RTFLUX data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE RTFLUX_PRINT(USER_RTFLUX)
IMPLICIT NONE
! Passed in
TYPE (RTFLUX_DATA) USER_RTFLUX  ! A user data variable to be defined by WRITEing in a RTFLUX file
! Local
DIF3D_Int I,J,K,G

100 FORMAT('[RTFLUX]...SORRY, BUT I MUST STOP')
101 FORMAT('[RTFLUX]',107('.'))
102 FORMAT('')
105 FORMAT('[RTFLUX]...THERE WAS A NON-FATAL ERROR THAT OCCURED IN (PRINT)',59('.'))

IF (.NOT. USER_RTFLUX%DEFINED) THEN  ! Data to write doesn't exist
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[RTFLUX]...USER DATA WAS UNDEFINED FOR PRINT <",A60,">",8("."))') USER_RTFLUX%FILENAME
   WRITE(MODULE_OUT,100)
   CALL Abort
ELSE
   WRITE(MODULE_OUT,'("[RTFLUX]...HNAME............................",55("."),A16)') USER_RTFLUX%HNAME
   WRITE(MODULE_OUT,'("[RTFLUX]...HUSE(1)..........................",55("."),A16)') USER_RTFLUX%HUSE(1)
   WRITE(MODULE_OUT,'("[RTFLUX]...HUSE(2)..........................",55("."),A16)') USER_RTFLUX%HUSE(2)
   WRITE(MODULE_OUT,'("[RTFLUX]...IVERS............................",55("."),I16)') USER_RTFLUX%IVERS
   WRITE(MODULE_OUT,'("[RTFLUX]...NUMBER OF DIMENSIONS.............",55("."),I16)') USER_RTFLUX%NDIM
   WRITE(MODULE_OUT,'("[RTFLUX]...NUMBER OF GROUPS.................",55("."),I16)') USER_RTFLUX%NGROUP
   WRITE(MODULE_OUT,'("[RTFLUX]...FIRST DIMENSION FINE INTERVALS...",55("."),I16)') USER_RTFLUX%NINTI
   WRITE(MODULE_OUT,'("[RTFLUX]...SECOND DIMENSION FINE INTERVALS..",55("."),I16)') USER_RTFLUX%NINTJ
   WRITE(MODULE_OUT,'("[RTFLUX]...THIRD DIMENSION FINE INTERVALS...",55("."),I16)') USER_RTFLUX%NINTK
   WRITE(MODULE_OUT,'("[RTFLUX]...EXPORTED OUTER ITERATION NUMBER..",55("."),I16)') USER_RTFLUX%ITER
   WRITE(MODULE_OUT,'("[RTFLUX]...EXPORTED K-EFFECTIVE.............",55("."),1PE16.9)') USER_RTFLUX%EFFK
   WRITE(MODULE_OUT,'("[RTFLUX]...NEUTRONICS POWER LEVEL (WATTS)...",55("."),1PE16.9)') USER_RTFLUX%POWER
   WRITE(MODULE_OUT,'("[RTFLUX]...DATA BLOCKING OF SECOND DIMENSION",55("."),I16)') USER_RTFLUX%NBLOK
   DO G = 1,USER_RTFLUX%NGROUP
      DO K = 1,USER_RTFLUX%NINTK
         WRITE(MODULE_OUT,'("[RTFLUX]",26("."),"REGULAR FLUX FOR AXIAL PLANE ",I3," GROUP ",26("."),I16)') K,G
         !WRITE(MODULE_OUT,'(("[RTFLUX]...J ->",3X, 5000(5X,I4,5X)))') (J,J=1,USER_RTFLUX%NINTJ)
         !DO I = 1,USER_RTFLUX%NINTI
         !   WRITE(MODULE_OUT,'(("[RTFLUX]...I=",I3,2X,5000(1PE13.5,1X)))') I,(USER_RTFLUX%FREG(I,J,K,G),J=1,USER_RTFLUX%NINTJ)
         !END DO
         ! This prints out the map like you want to see it
         WRITE(MODULE_OUT,'(("[RTFLUX]...J ->",3X, 5000(5X,I4,5X)))') (J,J=1,USER_RTFLUX%NINTJ)
         DO I = USER_RTFLUX%NINTI,1,-1
            WRITE(MODULE_OUT,'(("[RTFLUX]...I=",I3,2X,5000(1PE13.5,1X)))') I,(USER_RTFLUX%FREG(I,J,K,G),J=1,USER_RTFLUX%NINTJ)
         END DO
      END DO
   END DO
END IF

END SUBROUTINE RTFLUX_PRINT

!---------------------------------------------------------------------------------------------------------------------------------
!  RTFLUX_DEFINE(USER_RTFLUX,NINTI,NINTJ,NINTK,NGROUP)
!  Allocates and initializes the RTFLUX data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE RTFLUX_DEFINE(USER_RTFLUX,NINTI,NINTJ,NINTK,NGROUP)
IMPLICIT NONE
! Passed in
TYPE (RTFLUX_DATA) USER_RTFLUX    ! A user data variable to be defined by reading in a RTFLUX file
DIF3D_Int :: NINTI                       ! Number of increments in the I direction (X,R,etc...)
DIF3D_Int :: NINTJ                       ! Number of increments in the J direction (Y,THETA,etc...)
DIF3D_Int :: NINTK                       ! Number of increments in the K direction (Z,PHI,etc...)
DIF3D_Int :: NGROUP                      ! Number of groups
! Local
DIF3D_Int IOS

100 FORMAT('[RTFLUX]...SORRY, BUT I MUST STOP')
101 FORMAT('[RTFLUX]',107('.'))
102 FORMAT('')
105 FORMAT('[RTFLUX]...THERE WAS A FATAL ERROR THAT OCCURED IN (ALLOCATE)',54('.'))

IF (USER_RTFLUX%DEFINED) CALL RTFLUX_VOID(USER_RTFLUX)   ! If the variable is already defined, void it

#ifdef DIF3D_Debug
   WRITE(MODULE_OUT,'("[RTFLUX]...NUMBER OF ENERGY GROUPS........",57("."),I16)') NGROUP
   WRITE(MODULE_OUT,'("[RTFLUX]...FIRST  DIMENSION FINE INTERVALS",57("."),I16)') NINTI 
   WRITE(MODULE_OUT,'("[RTFLUX]...SECOND DIMENSION FINE INTERVALS",57("."),I16)') NINTJ 
   WRITE(MODULE_OUT,'("[RTFLUX]...THIRD  DIMENSION FINE INTERVALS",57("."),I16)') NINTK 
#endif

IF (NGROUP .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[RTFLUX]...INVALID NUMBER OF ENERGY GROUPS..................",39("."),I16)') NGROUP
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (NINTI .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[RTFLUX]...INVALID NUMBER OF FIRST  DIMENSION FINE INTERVALS",39("."),I16)') NINTI
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (NINTJ .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[RTFLUX]...INVALID NUMBER OF SECOND DIMENSION FINE INTERVALS",39("."),I16)') NINTJ
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (NINTK .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[RTFLUX]...INVALID NUMBER OF THIRD  DIMENSION FINE INTERVALS",39("."),I16)') NINTK
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

! Data for card type 0-4
ALLOCATE(USER_RTFLUX%FREG(NINTI,NINTJ,NINTK,NGROUP),STAT = IOS)
IF (IOS .NE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[RTFLUX]...FAILED TO ALLOCATE RTFLUX DATA",74("."))')
   WRITE(MODULE_OUT,'("[RTFLUX]...NUMBER OF ENERGY GROUPS........",57("."),I16)') NGROUP
   WRITE(MODULE_OUT,'("[RTFLUX]...FIRST  DIMENSION FINE INTERVALS",57("."),I16)') NINTI 
   WRITE(MODULE_OUT,'("[RTFLUX]...SECOND DIMENSION FINE INTERVALS",57("."),I16)') NINTJ 
   WRITE(MODULE_OUT,'("[RTFLUX]...THIRD  DIMENSION FINE INTERVALS",57("."),I16)') NINTK 
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

USER_RTFLUX%DEFINED  = .TRUE.
USER_RTFLUX%NGROUP   = NGROUP
USER_RTFLUX%NINTI    = NINTI
USER_RTFLUX%NINTJ    = NINTJ
USER_RTFLUX%NINTK    = NINTK
USER_RTFLUX%FREG     = 0.0d0

END SUBROUTINE RTFLUX_DEFINE

!---------------------------------------------------------------------------------------------------------------------------------
!  RTFLUX_VOID(USER_RTFLUX)
!  Deallocates the RTFLUX data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE RTFLUX_VOID(USER_RTFLUX)
IMPLICIT NONE
TYPE (RTFLUX_DATA) USER_RTFLUX  ! A user data variable to be defined by reading in a RTFLUX file
! local
DIF3D_Int IOS

100 FORMAT('[RTFLUX]...SORRY, BUT I MUST STOP')
101 FORMAT('[RTFLUX]',107('.'))
102 FORMAT('')
105 FORMAT('[RTFLUX]...THERE WAS A FATAL ERROR THAT OCCURED IN (VOIDTYPE)',54('.'))

IF (USER_RTFLUX%DEFINED) THEN
   DEALLOCATE(USER_RTFLUX%FREG,STAT = IOS)
   IF (IOS .NE. 0) THEN
      WRITE(MODULE_OUT,105)
      WRITE(MODULE_OUT,'("[RTFLUX]...FAILED TO DEALLOCATE RTFLUX DATA",72("."))')
      WRITE(MODULE_OUT,100)
      CALL Abort
   END IF
   USER_RTFLUX%DEFINED = .FALSE.
END IF

END SUBROUTINE RTFLUX_VOID

END MODULE RTFLUX_io
