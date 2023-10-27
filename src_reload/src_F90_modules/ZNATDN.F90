!---------------------------------------------------------------------------------------------------------------------------------
!  Using the historical format for the ZNATDN file, import the data from the binary file
!---------------------------------------------------------------------------------------------------------------------------------
!  ZNATDN_ASSIGNPRINTINFO     !  Sets the output unit and debugprint level settings for the module
!  ZNATDN_Copy                !  Duplicates a ZNATDN data structure
!  ZNATDN_DEFINE              !  Defines a ZNATDN_DATA type based upon user input
!  ZNATDN_VOID                !  Provides a path to deallocate the data
!  ZNATDN_PRINT               !  Prints a ZNATDN data structure
!  ZNATDN_MergeWithFactors    !  Merges OutZNATDN = F1*S1+F2*S2 where S1 and S2 are ZNATDN files and F1 and F2 are real scalars 
! ----------------------------
!  ZNATDN_ASSIGNPRINTINFO(OUTPUT_UNIT,DEBUGPRINT)
!  ZNATDN_Copy(Source_ZNATDN,Destination_ZNATDN)
!  ZNATDN_DEFINE(USER_ZNATDN,NNS,NTZSZ)
!  ZNATDN_VOID(USER_ZNATDN)
!  ZNATDN_PRINT(USER_ZNATDN)
!  ZNATDN_MergeWithFactors(S1ZNATDN,F1,S2ZNATDN,F2,OutZNATDN)
!---------------------------------------------------------------------------------------------------------------------------------
MODULE ZNATDN_io
#include "DIF3D_Types.h"
IMPLICIT NONE

TYPE ZNATDN_DATA
   DIF3D_Log  :: DEFINED  = .FALSE.
   DIF3D_FileNames       :: FILENAME = 'UNDEFINE'  ! Originating file name
   ! CARD TYPE 0
   DIF3D_Char        :: HNAME    = 'UNDEFINE'
   DIF3D_Char        :: HUSE(2)  = 'UNDEFINE'
   ZNATDN_Int :: IVERS    = 0
   ! CARD TYPE 1
   ZNATDN_Real :: TIME     = 0.0         ! REFERENCE REAL TIME, DAYS
   ZNATDN_Int :: NCY      = 0           ! REFERENCE CYCLE NUMBER
   ZNATDN_Int :: NTZSZ    = 0           ! NUMBER OF ZONES PLUS NUMBER OF SUBZONES
   ZNATDN_Int :: NNS      = 0           ! MAXIMUM NUMBER OF NUCLIDES IN ANY SET
   ZNATDN_Int :: NBLKAD   = 0           ! NUMBER OF BLOCKS OF ATOM DENSITY DATA 
   ! CARD TYPE 2
   ZNATDN_Real, DIMENSION(:,:), POINTER :: ADEN              ! NNS,NTZSZ
END TYPE

TYPE (ZNATDN_DATA), SAVE :: MODULE_ZNATDN ! A optional use "common block" data type for the module

DIF3D_Int            :: MODULE_OUT = 6               ! Module output unit
PRIVATE  MODULE_OUT          

CONTAINS

!---------------------------------------------------------------------------------------------------------------------------------
!  ZNATDN_ASSIGNPRINTINFO(OUTPUT_UNIT)
!  Sets the output unit and debugprint level settings for the module
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ZNATDN_ASSIGNPRINTINFO(MODULE_OUTPUT_UNIT)
DIF3D_Int MODULE_OUTPUT_UNIT    ! Reset the object output unit

MODULE_OUT = MODULE_OUTPUT_UNIT

END SUBROUTINE ZNATDN_ASSIGNPRINTINFO

!---------------------------------------------------------------------------------------------------------------------------------
!  ZNATDN_PRINT(USER_ZNATDN)
!  Prints a ZNATDN data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ZNATDN_PRINT(USER_ZNATDN)
IMPLICIT NONE
! Passed in
TYPE (ZNATDN_DATA) USER_ZNATDN  ! A user data variable to be defined by reading in a ZNATDN file
! Local
DIF3D_Int IOS,I,J,M,N

100 FORMAT('[ZNATDN]...Sorry, but I must stop')
101 FORMAT('[ZNATDN]',107('.'))
102 FORMAT('')
105 FORMAT('[ZNATDN]...THERE WAS A NON-FATAL ERROR THAT OCCURED IN (PRINT)',53('.'))

IF (.NOT. USER_ZNATDN%DEFINED) THEN  ! Data to write doesn't exist
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[ZNATDN]...USER DATA WAS UNDEFINED FOR PRINT <",A60,">",8("."))') USER_ZNATDN%FILENAME
   WRITE(MODULE_OUT,100)
ELSE
   WRITE(MODULE_OUT,'("[ZNATDN]...HNAME................................",51("."),A16)') USER_ZNATDN%HNAME
   WRITE(MODULE_OUT,'("[ZNATDN]...HUSE(1)..............................",51("."),A16)') USER_ZNATDN%HUSE(1)
   WRITE(MODULE_OUT,'("[ZNATDN]...HUSE(2)..............................",51("."),A16)') USER_ZNATDN%HUSE(2)
   WRITE(MODULE_OUT,'("[ZNATDN]...IVERS................................",51("."),I16)') USER_ZNATDN%IVERS
   WRITE(MODULE_OUT,'("[ZNATDN]...REFERENCE REAL TIME, DAYS............",51("."),1PE16.9)') USER_ZNATDN%TIME
   WRITE(MODULE_OUT,'("[ZNATDN]...REFERENCE CYCLE NUMBER...............",51("."),I16)') USER_ZNATDN%NCY
   WRITE(MODULE_OUT,'("[ZNATDN]...NUMBER OF ZONES PLUS SUBZONES........",51("."),I16)') USER_ZNATDN%NTZSZ
   WRITE(MODULE_OUT,'("[ZNATDN]...MAXIMUM NUMBER OF NUCLIDES IN ANY SET",51("."),I16)') USER_ZNATDN%NNS
   WRITE(MODULE_OUT,'("[ZNATDN]...NUMBER OF BLOCKS OF ATOM DENSITY DATA",51("."),I16)') USER_ZNATDN%NBLKAD
   DO J = 1,USER_ZNATDN%NTZSZ
      WRITE(MODULE_OUT,'("[ZNATDN]...ZONE ATOM DENSITIES FOR MIXTURE",57("."),I16)') J
      WRITE(MODULE_OUT,'(("[ZNATDN]...",4(I5," [",1PE13.6,"]",3X)))') (N,USER_ZNATDN%ADEN(N,J),N=1,USER_ZNATDN%NNS)
   END DO
END IF

END SUBROUTINE ZNATDN_PRINT

!---------------------------------------------------------------------------------------------------------------------------------
!  ZNATDN_Copy(Source_ZNATDN,Destination_ZNATDN)
!  Duplicates a ZNATDN data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ZNATDN_Copy(Source_ZNATDN,Destination_ZNATDN)
IMPLICIT NONE
! Passed in
TYPE (ZNATDN_DATA) Source_ZNATDN,Destination_ZNATDN
DIF3D_Int :: NNS           ! Number of nuclides in the data set
DIF3D_Int :: NTZSZ            ! Number of zones per subzone
! Local
DIF3D_Int I,J


CALL ZNATDN_Void(Destination_ZNATDN)
IF (Source_ZNATDN%Defined) THEN
   NNS = Source_ZNATDN%NNS
   NTZSZ  = Source_ZNATDN%NTZSZ
   CALL ZNATDN_DEFINE(Destination_ZNATDN,NNS,NTZSZ)
   Destination_ZNATDN%FILENAME = Source_ZNATDN%FILENAME
   Destination_ZNATDN%HNAME    = Source_ZNATDN%HNAME
   Destination_ZNATDN%HUSE(1)  = Source_ZNATDN%HUSE(1)
   Destination_ZNATDN%HUSE(2)  = Source_ZNATDN%HUSE(2)
   Destination_ZNATDN%IVERS    = Source_ZNATDN%IVERS
   Destination_ZNATDN%TIME     = Source_ZNATDN%TIME
   Destination_ZNATDN%NCY      = Source_ZNATDN%NCY
   Destination_ZNATDN%NBLKAD   = Source_ZNATDN%NBLKAD

   DO I = 1,Source_ZNATDN%NTZSZ
      DO J = 1,Source_ZNATDN%NNS
         Destination_ZNATDN%ADEN(J,I) = Source_ZNATDN%ADEN(J,I)
      END DO
   END DO
END IF

END SUBROUTINE ZNATDN_Copy

!---------------------------------------------------------------------------------------------------------------------------------
!  ZNATDN_DEFINE(USER_ZNATDN,NNS,NTZSZ)
!  Allocates and initializes the ZNATDN data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ZNATDN_DEFINE(USER_ZNATDN,NNS,NTZSZ)
IMPLICIT NONE
! Passed in
TYPE (ZNATDN_DATA) USER_ZNATDN  ! A user data variable to be defined by reading in a ZNATDN file
DIF3D_Int :: NNS           ! Number of nuclides in the data set
DIF3D_Int :: NTZSZ            ! Number of zones per subzone
! Local
DIF3D_Int IOS

100 FORMAT('[ZNATDN]...Sorry, but I must stop')
101 FORMAT('[ZNATDN]',107('.'))
102 FORMAT('')
105 FORMAT('[ZNATDN]...There was a fatal error that occured in (ALLOCATE)',54('.'))

IF (USER_ZNATDN%DEFINED) CALL ZNATDN_VOID(USER_ZNATDN)   ! If the variable is already defined, void it

#ifdef DIF3D_Debug
   WRITE(MODULE_OUT,'("[ZNATDN]...MAXIMUM NUMBER OF NUCLIDES IN ANY SET",51("."),I16)') NNS 
   WRITE(MODULE_OUT,'("[ZNATDN]...NUMBER OF ZONES PLUS SUBZONES",59("."),I16)') NTZSZ
#endif

IF (NNS .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[ZNATDN]...INVALID NUMBER OF NUCLIDES IN A SET",53("."),I16)') NNS
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (NTZSZ .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[ZNATDN]...INVALID NUMBER OF ZONES PLUS SUBZONES",51("."),I16)') NTZSZ
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

! Data for card type 0-4
ALLOCATE(USER_ZNATDN%ADEN(NNS,NTZSZ),STAT = IOS)
IF (IOS .NE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[ZNATDN]...FAILED TO ALLOCATE ZNATDN DATA",74("."))')
   WRITE(MODULE_OUT,'("[ZNATDN]...MAXIMUM NUMBER OF NUCLIDES IN ANY SET",51("."),I16)') NNS 
   WRITE(MODULE_OUT,'("[ZNATDN]...NUMBER OF ZONES PLUS SUBZONES",59("."),I16)') NTZSZ
   WRITE(MODULE_OUT,100)
   CALL Basic_Abort
END IF

USER_ZNATDN%DEFINED  = .TRUE.
USER_ZNATDN%NNS      = NNS
USER_ZNATDN%NTZSZ    = NTZSZ
USER_ZNATDN%ADEN     = 0.0

END SUBROUTINE ZNATDN_DEFINE

!---------------------------------------------------------------------------------------------------------------------------------
!  ZNATDN_VOID(USER_ZNATDN)
!  Deallocates the data
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ZNATDN_VOID(USER_ZNATDN)
IMPLICIT NONE
TYPE (ZNATDN_DATA) USER_ZNATDN  ! A user data variable to be defined by reading in a ZNATDN file
! local
DIF3D_Int IOS

100 FORMAT('[ZNATDN]...Sorry, but I must stop')
101 FORMAT('[ZNATDN]',107('.'))
102 FORMAT('')
105 FORMAT('[ZNATDN]...There was a fatal error that occured in (VOIDTYPE)',54('.'))

IF (USER_ZNATDN%DEFINED) THEN
   DEALLOCATE(USER_ZNATDN%ADEN,STAT = IOS)
   IF (IOS .NE. 0) THEN
      WRITE(MODULE_OUT,105)
      WRITE(MODULE_OUT,'("[ZNATDN]...FAILED TO DEALLOCATE ZNATDN DATA",72("."))')
      WRITE(MODULE_OUT,100)
      CALL Basic_Abort
   END IF
   USER_ZNATDN%DEFINED = .FALSE.
END IF

END SUBROUTINE ZNATDN_VOID

!---------------------------------------------------------------------------------------------------------------------------------
!  ZNATDN_MergeWithFactors(S1ZNATDN,F1,S2ZNATDN,F2,OutZNATDN)
!  Merges OutZNATDN = F1*S1+F2*S2 where S1 and S2 are ZNATDN files and F1 and F2 are real scalars 
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ZNATDN_MergeWithFactors(S1,F1,S2,F2,OutZNATDN)
IMPLICIT NONE
! Passed in
TYPE (ZNATDN_DATA) S1,S2,OutZNATDN
DIF3D_Real F1,F2
! Local
DIF3D_Int I,J

100 FORMAT('[ZNATDN]...Sorry, but I must stop')
101 FORMAT('[ZNATDN]',107('.'))
102 FORMAT('')
105 FORMAT('[ZNATDN]...There was a fatal error that occured in (MergeWithFactors)',54('.'))

CALL ZNATDN_Void(OutZNATDN)
IF ((S1%Defined) .AND. (S2%Defined)) THEN
   IF (S1%NTZSZ .NE. S2%NTZSZ) THEN
      WRITE(MODULE_OUT,105)
      WRITE(MODULE_OUT,'("[ZNATDN]...S1 and S2 have different NTZSZ ")')
      WRITE(MODULE_OUT,100)
      CALL Basic_Abort
   ELSE IF (S1%NNS .NE. S2%NNS) THEN
      WRITE(MODULE_OUT,105)
      WRITE(MODULE_OUT,'("[ZNATDN]...S1 and S2 have different NNS ")')
      WRITE(MODULE_OUT,100)
      CALL Basic_Abort
   ELSE
      CALL ZNATDN_Copy(S1,OutZNATDN)
      DO I = 1,S1%NTZSZ
         DO J = 1,S1%NNS
            OutZNATDN%ADEN(J,I) = F1*S1%ADEN(J,I) + F2*S2%ADEN(J,I)
         END DO
      END DO
   END IF
END IF

END SUBROUTINE ZNATDN_MergeWithFactors

END MODULE ZNATDN_io
