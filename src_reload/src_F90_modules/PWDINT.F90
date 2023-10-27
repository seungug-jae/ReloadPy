!---------------------------------------------------------------------------------------------------------------------------------
!  Using the historical format for the PWDINT file, import/export data from the binary file 
!---------------------------------------------------------------------------------------------------------------------------------
!  PWDINT_ASSIGNPRINTINFO   !  Sets the output unit and debugprint level settings for the module
!  PWDINT_TotalPower        !  Computes the total power of the system
!  PWDINT_DEFINE            !  Defines a PWDINT_DATA type based upon user input
!  PWDINT_VOID              !  Provides a path to deallocate the data
!  PWDINT_PRINT             !  Prints the PWDINT data structure
!---------------------------
!  PWDINT_ASSIGNPRINTINFO(OUTPUT_UNIT)
!  PWDINT_TotalPower(PWDINT,TotalPower,MeshVolumes,MeshFluxFactors)
!  PWDINT_DEFINE(USER_PWDINT,USER_NINTI,USER_NINTJ,USER_NINTK)
!  PWDINT_VOID(USER_PWDINT)
!  PWDINT_PRINT(USER_PWDINT)
!---------------------------------------------------------------------------------------------------------------------------------
MODULE PWDINT_io
#include "DIF3D_Types.h"
IMPLICIT NONE


TYPE PWDINT_DATA
   DIF3D_Log       :: DEFINED=.FALSE.
   DIF3D_FileNames       :: FILENAME = 'UNDEFINE'  ! Originating file name
   ! CARD TYPE 0
   DIF3D_Char        :: HNAME    = 'UNDEFINE'
   DIF3D_Char        :: HUSE(2)  = 'UNDEFINE'
   PWDINT_Int :: IVERS    = 0
   ! CARD TYPE 1
   PWDINT_Real :: TIME  = 0.0            ! REFERENCE REAL TIME, DAYS                 
   PWDINT_Real :: POWER = 0.0            ! POWER LEVEL FOR ACTUAL NEUTRONICS PROBLEM, WATTS THERMAL
   PWDINT_Real :: VOL   = 0.0            ! VOLUME (CC) OVER WHICH POWER WAS DETERMINED 
   PWDINT_Int :: NINTI = 0              ! NUMBER OF FIRST DIMENSION FINE INTERVALS  
   PWDINT_Int :: NINTJ = 0              ! NUMBER OF SECOND DIMENSION FINE INTERVALS 
   PWDINT_Int :: NINTK = 0              ! NUMBER OF THIRD DIMENSION FINE INTERVALS  
   PWDINT_Int :: NCY   = 0              ! REFERENCE COUNT (CYCLE NUMBER)            
   PWDINT_Int :: NBLOK = 0              ! DATA BLOCKING FACTOR. THE SECOND DIMENSION VARIABLE IS BLOCKED INTO NBLOK BLOCKS.
   ! CARD TYPE 2
   PWDINT_Real, DIMENSION(:,:,:), POINTER :: PWR  ! POWER DENSITY BY INTERVAL, WATTS/CC (NINTI,NINTJ,NINTK)
END TYPE

TYPE (PWDINT_DATA), SAVE :: MODULE_PWDINT          ! A optional use "common block" data type for the module

DIF3D_Int :: MODULE_OUT = 6               ! Module output unit
PRIVATE  MODULE_OUT

CONTAINS

!---------------------------------------------------------------------------------------------------------------------------------
!  PWDINT_ASSIGNPRINTINFO(OUTPUT_UNIT)
!  Sets the output unit and debugprint level settings for the module
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE PWDINT_ASSIGNPRINTINFO(MODULE_OUTPUT_UNIT)
DIF3D_Int MODULE_OUTPUT_UNIT    ! Reset the object output unit

MODULE_OUT = MODULE_OUTPUT_UNIT

END SUBROUTINE PWDINT_ASSIGNPRINTINFO

!---------------------------------------------------------------------------------------------------------------------------------
!  PWDINT_TotalPower(PWDINT,TotalPower,MeshVolumes,MeshFluxFactors)
!  Computes the total power of the system
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE PWDINT_TotalPower(PWDINT,TotalPower,MeshVolumes,MeshFluxFactors)
TYPE (PWDINT_DATA) :: PWDINT
DIF3D_Real         :: TotalPower
DIF3D_Real         :: MeshVolumes(PWDINT%NINTI,PWDINT%NINTJ,PWDINT%NINTK)  ! GEODST_MeshVolumes
DIF3D_Real         :: MeshFluxFactors(PWDINT%NINTI,PWDINT%NINTJ)           ! GEODST_NODAL_factors
! Local
DIF3D_Int I,J,K

TotalPower = 0.0D0
DO K = 1,PWDINT%NINTK
   DO J = 1,PWDINT%NINTJ
      DO I = 1,PWDINT%NINTI
         TotalPower = TotalPower + PWDINT%PWR(I,J,K) * MeshVolumes(I,J,K) * MeshFluxFactors(I,J)
      END DO
   END DO
END DO

END SUBROUTINE PWDINT_TotalPower

!---------------------------------------------------------------------------------------------------------------------------------
!  PWDINT_PRINT(USER_PWDINT)
!  Prints the PWDINT data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE PWDINT_PRINT(USER_PWDINT)
IMPLICIT NONE
! Passed in
TYPE (PWDINT_DATA) USER_PWDINT  ! A user data variable to be defined by WRITEing in a PWDINT file
! Local
DIF3D_Int I,J,K

100 FORMAT('[PWDINT]...SORRY, BUT I MUST STOP')
101 FORMAT('[PWDINT]',107('.'))
102 FORMAT('')
105 FORMAT('[PWDINT]...THERE WAS A NON-FATAL ERROR THAT OCCURED IN (PRINT)',59('.'))

IF (.NOT. USER_PWDINT%DEFINED) THEN  ! Data to write doesn't exist
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[PWDINT]...USER DATA WAS UNDEFINED FOR PRINT <",A60,">",8("."))') USER_PWDINT%FILENAME
   WRITE(MODULE_OUT,100)
   CALL Abort
ELSE
   WRITE(MODULE_OUT,'("[PWDINT]...HNAME............................",55("."),A16)') USER_PWDINT%HNAME
   WRITE(MODULE_OUT,'("[PWDINT]...HUSE(1)..........................",55("."),A16)') USER_PWDINT%HUSE(1)
   WRITE(MODULE_OUT,'("[PWDINT]...HUSE(2)..........................",55("."),A16)') USER_PWDINT%HUSE(2)
   WRITE(MODULE_OUT,'("[PWDINT]...IVERS............................",55("."),I16)') USER_PWDINT%IVERS
   WRITE(MODULE_OUT,'("[PWDINT]...REFERENCE REAL TIME, DAYS........",55("."),1PE16.9)') USER_PWDINT%TIME
   WRITE(MODULE_OUT,'("[PWDINT]...NEUTRONICS POWER LEVEL (WATTS)...",55("."),1PE16.9)') USER_PWDINT%POWER
   WRITE(MODULE_OUT,'("[PWDINT]...VOLUME (CC) OF INTERVAL..........",55("."),1PE16.9)') USER_PWDINT%VOL
   WRITE(MODULE_OUT,'("[PWDINT]...FIRST DIMENSION FINE INTERVALS...",55("."),I16)') USER_PWDINT%NINTI
   WRITE(MODULE_OUT,'("[PWDINT]...SECOND DIMENSION FINE INTERVALS..",55("."),I16)') USER_PWDINT%NINTJ
   WRITE(MODULE_OUT,'("[PWDINT]...THIRD DIMENSION FINE INTERVALS...",55("."),I16)') USER_PWDINT%NINTK
   WRITE(MODULE_OUT,'("[PWDINT]...REFERENCE COUNT (CYCLE NUMBER)...",55("."),I16)') USER_PWDINT%NCY
   WRITE(MODULE_OUT,'("[PWDINT]...DATA BLOCKING OF SECOND DIMENSION",55("."),I16)') USER_PWDINT%NBLOK
   DO K = 1,USER_PWDINT%NINTK
      WRITE(MODULE_OUT,'("[PWDINT]",31("."),"POWER DENSITY FOR AXIAL PLANE",31("."),I16)') K
         WRITE(MODULE_OUT,'(("[PWDINT]...J ->",3X, 5000(5X,I4,5X)))') (J,J=1,USER_PWDINT%NINTJ)
      !DO I = 1,USER_PWDINT%NINTI
      !   WRITE(MODULE_OUT,'(("[PWDINT]...I=",I3,2X,5000(1PE13.5,1X)))') I,(USER_PWDINT%PWR(I,J,K),J=1,USER_PWDINT%NINTJ)
      !END DO
      DO I = USER_PWDINT%NINTI,1,-1
         WRITE(MODULE_OUT,'(("[PWDINT]...I=",I3,2X,5000(1PE13.5,1X)))') I,(USER_PWDINT%PWR(I,J,K),J=1,USER_PWDINT%NINTJ)
      END DO
   END DO
END IF

END SUBROUTINE PWDINT_PRINT

!---------------------------------------------------------------------------------------------------------------------------------
!  PWDINT_DEFINE(USER_PWDINT,NINTI,NINTJ,NINTK)
!  Allocates and initializes the PWDINT data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE PWDINT_DEFINE(USER_PWDINT,NINTI,NINTJ,NINTK)
IMPLICIT NONE
! Passed in
TYPE (PWDINT_DATA) USER_PWDINT  ! A user data variable to be defined by reading in a PWDINT file
DIF3D_Int :: NINTI                      ! Number of increments in the I direction (X,R,etc...)
DIF3D_Int :: NINTJ                      ! Number of increments in the J direction (Y,THETA,etc...)
DIF3D_Int :: NINTK                      ! Number of increments in the K direction (Z,PHI,etc...)
! Local
DIF3D_Int IOS

100 FORMAT('[PWDINT]...SORRY, BUT I MUST STOP')
101 FORMAT('[PWDINT]',107('.'))
102 FORMAT('')
105 FORMAT('[PWDINT]...THERE WAS A FATAL ERROR THAT OCCURED IN (ALLOCATE)',54('.'))

IF (USER_PWDINT%DEFINED) CALL PWDINT_VOID(USER_PWDINT)   ! If the variable is already defined, void it

#ifdef DIF3D_Debug
   WRITE(MODULE_OUT,'("[PWDINT]...FIRST  DIMENSION FINE INTERVALS",57("."),I16)') NINTI 
   WRITE(MODULE_OUT,'("[PWDINT]...SECOND DIMENSION FINE INTERVALS",57("."),I16)') NINTJ 
   WRITE(MODULE_OUT,'("[PWDINT]...THIRD  DIMENSION FINE INTERVALS",57("."),I16)') NINTK 
#endif

IF (NINTI .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[PWDINT]...INVALID NUMBER OF FIRST  DIMENSION FINE INTERVALS",39("."),I16)') NINTI
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (NINTJ .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[PWDINT]...INVALID NUMBER OF SECOND DIMENSION FINE INTERVALS",39("."),I16)') NINTJ
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (NINTK .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[PWDINT]...INVALID NUMBER OF THIRD  DIMENSION FINE INTERVALS",39("."),I16)') NINTK
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

! Data for card type 0-4
ALLOCATE(USER_PWDINT%PWR(NINTI,NINTJ,NINTK),STAT = IOS)
IF (IOS .NE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[PWDINT]...FAILED TO ALLOCATE PWDINT DATA",74("."))')
   WRITE(MODULE_OUT,'("[PWDINT]...FIRST  DIMENSION FINE INTERVALS",57("."),I16)') NINTI 
   WRITE(MODULE_OUT,'("[PWDINT]...SECOND DIMENSION FINE INTERVALS",57("."),I16)') NINTJ 
   WRITE(MODULE_OUT,'("[PWDINT]...THIRD  DIMENSION FINE INTERVALS",57("."),I16)') NINTK 
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

USER_PWDINT%DEFINED  = .TRUE.
USER_PWDINT%NINTI    = NINTI
USER_PWDINT%NINTJ    = NINTJ
USER_PWDINT%NINTK    = NINTK
USER_PWDINT%PWR      = 0.0

END SUBROUTINE PWDINT_DEFINE

!---------------------------------------------------------------------------------------------------------------------------------
!  PWDINT_VOID(USER_PWDINT)
!  Deallocates the PWDINT data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE PWDINT_VOID(USER_PWDINT)
IMPLICIT NONE
TYPE (PWDINT_DATA) USER_PWDINT  ! A user data variable to be defined by reading in a PWDINT file
! local
DIF3D_Int IOS

100 FORMAT('[PWDINT]...SORRY, BUT I MUST STOP')
101 FORMAT('[PWDINT]',107('.'))
102 FORMAT('')
105 FORMAT('[PWDINT]...THERE WAS A FATAL ERROR THAT OCCURED IN (VOIDTYPE)',54('.'))

IF (USER_PWDINT%DEFINED) THEN
   DEALLOCATE(USER_PWDINT%PWR,STAT = IOS)
   IF (IOS .NE. 0) THEN
      WRITE(MODULE_OUT,105)
      WRITE(MODULE_OUT,'("[PWDINT]...FAILED TO DEALLOCATE PWDINT DATA",72("."))')
      WRITE(MODULE_OUT,100)
      CALL Abort
   END IF
   USER_PWDINT%DEFINED = .FALSE.
END IF

END SUBROUTINE PWDINT_VOID

END MODULE PWDINT_io
