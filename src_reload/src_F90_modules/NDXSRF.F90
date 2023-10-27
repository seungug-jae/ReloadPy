!---------------------------------------------------------------------------------------------------------------------------------
!  Using the historical format for the NDXSRF file, import the data from the binary file 
!---------------------------------------------------------------------------------------------------------------------------------
!  NDXSRF_ASSIGNPRINTINFO  !  Sets the output unit and debugprint level settings for the module
!  NDXSRF_PRINT            !  Prints the NDXSRF data structure
!  NDXSRF_Copy             !  Duplicates the NDXSRF data structure
!  NDXSRF_DEFINE           !  Allocates and initializes the NDXSRF data structure
!  NDXSRF_VOID             !  Deallocates the data
!--------------------------
!  NDXSRF_ASSIGNPRINTINFO(OUTPUT_UNIT)
!  NDXSRF_PRINT(USER_NDXSRF)
!  NDXSRF_Copy(Source_NDXSRF,Destination_NDXSRF)
!  NDXSRF_DEFINE(USER_NDXSRF,NON,NSN,NNS,NAN,NZONE,NSZ)
!  NDXSRF_VOID(USER_NDXSRF)
!---------------------------------------------------------------------------------------------------------------------------------
MODULE NDXSRF_io 
#include "DIF3D_Types.h"
IMPLICIT NONE

TYPE NDXSRF_DATA
   DIF3D_Log       :: DEFINED  = .FALSE.
   DIF3D_FileNames :: FILENAME = 'UNDEFINE'  ! Originating file name
   ! CARD TYPE 0
   DIF3D_Char        :: HNAME    = 'UNDEFINE'
   DIF3D_Char        :: HUSE(2)  = 'UNDEFINE'
   NDXSRF_Int :: IVERS    = 0
   ! CARD TYPE 1
   NDXSRF_Int :: NON      = 0                      ! NUMBER OF NUCLIDES IN CROSS SECTION DATA
   NDXSRF_Int :: NSN      = 0                      ! NUMBER OF NUCLIDE SETS IDENTIFIED       
   NDXSRF_Int :: NNS      = 0                      ! MAXIMUM NUMBER OF NUCLIDES IN ANY SET   
   NDXSRF_Int :: NAN      = 0                      ! NUMBER OF DIFFERENT NUCLIDES IN DATA    
   NDXSRF_Int :: NZONE    = 0                      ! NUMBER OF ZONES                         
   NDXSRF_Int :: NSZ      = 0                      ! NUMBER OF SUBZONES (SUBASSEMBLIES)      
   ! CARD TYPE 2
   DIF3D_Char,       DIMENSION(:),  POINTER :: HNNAME ! NON       HNNAME(N)   UNIQUE REFERENCE NUCLIDE NAME, IN LIBRARY ORDER (A6) ALPHANUMERIC
   DIF3D_Char,       DIMENSION(:),  POINTER :: HANAME ! NON       HANAME(N)   ABSOLUTE NUCLIDE REFERENCE, IN LIBRARY ORDER (A6) ALPHANUMERIC
   NDXSRF_Real,DIMENSION(:),  POINTER :: WPF    ! NON       WPF(N)      RESERVED                                       
   NDXSRF_Real,DIMENSION(:),  POINTER :: ATWT   ! NAN       ATWT(J)     ATOMIC WEIGHT                                  
   NDXSRF_Int, DIMENSION(:),  POINTER :: NCLN   ! NON       NCLN(N)     CLASSIFICATION  1- FISSILE            2- FERTILE      3- OTHER ACTINIDE
                                                           !                                       4- FISSION PRODUCT    5- STRUCTURAL   6- COOLANT
                                                           !                                       7- CONTROL ROD        7- UNDEFINED
   NDXSRF_Int,DIMENSION(:,:),POINTER :: NDXS   ! 4,NSN     NDXS(K,L)   REFERENCE DATA FOR SET L   K = 1, # OF NUCLIDES IN SET   K = 2, RESERVED
                                                           !                                                  K = 3, RESERVED         K = 4, RESERVED
   NDXSRF_Int,DIMENSION(:,:),POINTER :: NOS    ! NNS,NSN   NOS(I,L)    ORDER # OF NUCLIDE IN CROSS SECTION DATA (IN HNNAME LIST) OF NUCLIDE I IN SET L
   NDXSRF_Int,DIMENSION(:,:),POINTER :: NOR    ! NON,NSN   NOR(N,L)    ORDER NUMBER OF NUCLIDE IN SET L GIVEN ORDER NUMBER N IN CROSS SECTION DATA
   ! CARD TYPE 3
   !    NOTE THAT TO CALCULATE MACROSCOPIC CROSS SECTIONS FOR A ZONE, IT IS NECESSARY TO CONSIDER THE CONCENTRATION OF EACH NUCLIDE    
   !    IN THE PRIMARY SET ASSIGNMENT (UNLESS A ZERO IN NSPA INDICATES THERE ARE NONE) TIMES THE VOLUME FRACTION, AND THE CONCENTRATION 
   !    OF EACH NUCLIDE IN EACH SUBZONE ASSIGNED TO THE ZONE TIMES THE RATIO OF THE SUBZONE VOLUME TO THE ZONE VOLUME.                  
   NDXSRF_Real,DIMENSION(:),POINTER :: VOLZ     ! NZONE    VOLZ(N) VOLUMES OF ZONES, CC                           
   NDXSRF_Real,DIMENSION(:),POINTER :: VFPA     ! NZONE    VFPA(N) VOLUME FRACTIONS FOR PRIMARY ZONE ASSIGNMENTS  
   NDXSRF_Real,DIMENSION(:),POINTER :: VLSA     ! NSZ      VLSA(M) VOLUMES OF SUBZONES                            
   NDXSRF_Int,DIMENSION(:),POINTER :: NSPA     ! NZONE    NSPA(N) NUCLIDE SET REFERENCE, PRIMARY ZONE ASSIGNMENT (MAY BE ZERO ONLY IF THERE ARE SUBZONES)
   NDXSRF_Int,DIMENSION(:),POINTER :: NSSA     ! NSZ      NSSA(M) NUCLIDE SET REFERENCE ASSIGNMENT TO SUBZONES   
   NDXSRF_Int,DIMENSION(:),POINTER :: NZSZ     ! NSZ      NZSZ(M) ZONE CONTAINING SUBZONE                        
   ! Understanding the above pile of shit which is this way because of REBUS
   ! Composition (COMPXS) = zone
   ! zone = sum of all isotopes and all isotopes in all "sub-zones" assigned to that zone via NZSZ
   ! For zone level isotopes, we just use the Some_XS * atom_density(ZNATDN) * VFPA
   ! For the sub-zone level isotopes, we use Some_XS * atom_density(ZNATDN) * VLSA/VOLZ
END TYPE

TYPE (NDXSRF_DATA), SAVE :: MODULE_NDXSRF ! A optional use "common block" data type for the module

DIF3D_Int            :: MODULE_OUT = 6               ! Module output unit

PRIVATE  MODULE_OUT                                  ! Variables

CONTAINS

!---------------------------------------------------------------------------------------------------------------------------------
!  NDXSRF_ASSIGNPRINTINFO(OUTPUT_UNIT)
!  Sets the output unit and debugprint level settings for the module
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE NDXSRF_ASSIGNPRINTINFO(MODULE_OUTPUT_UNIT)
DIF3D_Int MODULE_OUTPUT_UNIT    ! Reset the object output unit

MODULE_OUT = MODULE_OUTPUT_UNIT

END SUBROUTINE NDXSRF_ASSIGNPRINTINFO

!---------------------------------------------------------------------------------------------------------------------------------
!  NDXSRF_PRINT(USER_NDXSRF)
!  Prints the NDXSRF data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE NDXSRF_PRINT(USER_NDXSRF)
IMPLICIT NONE
! Passed in
TYPE (NDXSRF_DATA) USER_NDXSRF  ! A user data variable to be defined by WRITEing in a NDXSRF file
! Local
DIF3D_Int IOS,I,J,K,L,M,N

100 FORMAT('[NDXSRF]...SORRY, BUT I MUST STOP')
101 FORMAT('[NDXSRF]',107('.'))
102 FORMAT('')
105 FORMAT('[NDXSRF]...THERE WAS A NON-FATAL ERROR THAT OCCURED IN (PRINT)',53('.'))

IF (.NOT. USER_NDXSRF%DEFINED) THEN  ! Data to write doesn't exist
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[NDXSRF]...USER DATA WAS UNDEFINED FOR PRINT <",A60,">",8("."))') USER_NDXSRF%FILENAME
   WRITE(MODULE_OUT,100)
   CALL Abort
ELSE
   WRITE(MODULE_OUT,'("[NDXSRF]...HNAME...................................",48("."),A16)') USER_NDXSRF%HNAME
   WRITE(MODULE_OUT,'("[NDXSRF]...HUSE(1).................................",48("."),A16)') USER_NDXSRF%HUSE(1)
   WRITE(MODULE_OUT,'("[NDXSRF]...HUSE(2).................................",48("."),A16)') USER_NDXSRF%HUSE(2)
   WRITE(MODULE_OUT,'("[NDXSRF]...IVERS...................................",48("."),I16)') USER_NDXSRF%IVERS
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF NUCLIDES IN CROSS SECTION DATA",48("."),I16)') USER_NDXSRF%NON
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF NUCLIDE SETS IDENTIFIED.......",48("."),I16)') USER_NDXSRF%NSN
   WRITE(MODULE_OUT,'("[NDXSRF]...MAXIMUM NUMBER OF NUCLIDES IN ANY SET...",48("."),I16)') USER_NDXSRF%NNS
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF DIFFERENT NUCLIDES IN DATA....",48("."),I16)') USER_NDXSRF%NAN
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF ZONES.........................",48("."),I16)') USER_NDXSRF%NZONE
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF SUBZONES (SUBASSEMBLIES)......",48("."),I16)') USER_NDXSRF%NSZ  
   IF (USER_NDXSRF%NAN .EQ. USER_NDXSRF%NON) THEN
      DO I = 1,USER_NDXSRF%NON
         WRITE(MODULE_OUT,301)I,USER_NDXSRF%HNNAME(I),USER_NDXSRF%HANAME(I),USER_NDXSRF%WPF(I),USER_NDXSRF%NCLN(I),USER_NDXSRF%ATWT(I)
         301 FORMAT('[NDXSRF]...NUCLIDE [',I5,',',A8,',',A8,'] RESERVED [',1PE13.6,'] CLASSIFICATION [',I2,']',' ATOM WEIGHT ',0PF16.9)
      END DO
   ELSE
      DO I = 1,USER_NDXSRF%NON
         WRITE(MODULE_OUT,300)I,USER_NDXSRF%HNNAME(I),USER_NDXSRF%HANAME(I),USER_NDXSRF%WPF(I),USER_NDXSRF%NCLN(I)
         300 FORMAT('[NDXSRF]...NUCLIDE [',I5,',',A8,',',A8,'] RESERVED [',1PE13.6,'] CLASSIFICATION [',I2,']',26('.'))
      END DO
      DO I = 1,USER_NDXSRF%NAN
         WRITE(MODULE_OUT,'("[NDXSRF]...ATOM WEIGHT OF UNIQUE ISOTOPE ",I4,54("."),F16.9)') I,USER_NDXSRF%ATWT(I)
      END DO
   END IF
   DO J = 1,USER_NDXSRF%NSN
      WRITE(MODULE_OUT,'("[NDXSRF]...THERE ARE ",I5," NUCLIDES (NDXS) IN SET",50("."),I16)') USER_NDXSRF%NDXS(1,J),J
      !WRITE(MODULE_OUT,'("[NDXSRF]...NOS[I,#] ORDER # OF NUCLIDE I IN CROSS SECTION DATA (HNNAME) FOR THIS SET",34("."))')
      !WRITE(MODULE_OUT,'("[NDXSRF]...",8("[",I4,",",I4,"]",1X))') (I,USER_NDXSRF%NOS(I,J),I=1,USER_NDXSRF%NDXS(1,J))
      !WRITE(MODULE_OUT,'("[NDXSRF]...NOR[I,#] ORDER # IN THIS SET OF NUCLIDE I FROM CROSS SECTION DATA (HNNAME)")')
      !WRITE(MODULE_OUT,'("[NDXSRF]...",8("[",I4,",",I4,"]",1X))') (I,USER_NDXSRF%NOR(I,J),I=1,USER_NDXSRF%NON)
      DO I = 1,USER_NDXSRF%NNS ! USER_NDXSRF%NDXS(1,J) ! Will be less than USER_NDXSRF%NON
         WRITE(MODULE_OUT,'("[NDXSRF]...Index ",I6," ISOTXS ID ",I6," SET INDEX ID ",I6)') &
            I,USER_NDXSRF%NOS(I,J),USER_NDXSRF%NOR(I,J)
      END DO
      DO I = USER_NDXSRF%NNS+1,USER_NDXSRF%NON !USER_NDXSRF%NDXS(1,J)+1,USER_NDXSRF%NON
         WRITE(MODULE_OUT,'("[NDXSRF]...Index ",I6," --------- ",6X," SET INDEX ID ",I6)') &
            I,USER_NDXSRF%NOR(I,J)
      END DO
   END DO
   DO I = 1,USER_NDXSRF%NZONE
      WRITE(MODULE_OUT,350)I,USER_NDXSRF%VOLZ(I),USER_NDXSRF%VFPA(I),USER_NDXSRF%NSPA(I)
      350 FORMAT('[NDXSRF]...ZONE ',I5,' OF VOLUME [',1PE13.6,',',1PE13.6,']',' HAS A PRIMARY ZONE ASSIGNMENT OF ',4('.'),I16)
   END DO
   DO I = 1,USER_NDXSRF%NSZ
      WRITE(MODULE_OUT,400)I,USER_NDXSRF%VLSA(I),USER_NDXSRF%NSSA(I),USER_NDXSRF%NZSZ(I)
      400 FORMAT('[NDXSRF]...SUBZONE ',I5,' OF VOLUME ',1PE13.6,' IS IN SET ',I5,' AND IS PART OF ZONE',15('.'),I16)
   END DO
END IF

END SUBROUTINE NDXSRF_PRINT

!---------------------------------------------------------------------------------------------------------------------------------
!  NDXSRF_Copy(Source_NDXSRF,Destination_NDXSRF)
!  Duplicates the NDXSRF data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE NDXSRF_Copy(Source_NDXSRF,Destination_NDXSRF)
IMPLICIT NONE
! Passed in
TYPE (NDXSRF_DATA) Source_NDXSRF,Destination_NDXSRF 
! Local
DIF3D_Int :: NON                       ! Number of nuclides in the cross section data
DIF3D_Int :: NSN                       ! Number of nuclide sets (ususually 1)
DIF3D_Int :: NNS                       ! Maximum number of nuclides in a set
DIF3D_Int :: NAN                       ! Number of unique nuclides in the cross section data
DIF3D_Int :: NZONE                     ! Number of zones
DIF3D_Int :: NSZ                       ! Number of subzones
! Local
DIF3D_Int I,J

CALL NDXSRF_Void(Destination_NDXSRF)

IF (Source_NDXSRF%Defined) THEN
   NON   = Source_NDXSRF%NON  
   NSN   = Source_NDXSRF%NSN  
   NNS   = Source_NDXSRF%NNS  
   NAN   = Source_NDXSRF%NAN  
   NZONE = Source_NDXSRF%NZONE
   NSZ   = Source_NDXSRF%NSZ
   CALL NDXSRF_DEFINE(Destination_NDXSRF,NON,NSN,NNS,NAN,NZONE,NSZ) 
   Destination_NDXSRF%FILENAME = Source_NDXSRF%FILENAME
   Destination_NDXSRF%HNAME    = Source_NDXSRF%HNAME
   Destination_NDXSRF%HUSE(1)  = Source_NDXSRF%HUSE(1)
   Destination_NDXSRF%HUSE(2)  = Source_NDXSRF%HUSE(2)
   Destination_NDXSRF%IVERS    = Source_NDXSRF%IVERS
   DO I = 1,NON
      Destination_NDXSRF%HNNAME(I) = Source_NDXSRF%HNNAME(I)
      Destination_NDXSRF%HANAME(I) = Source_NDXSRF%HANAME(I)
      Destination_NDXSRF%WPF(I)    = Source_NDXSRF%WPF(I)
      Destination_NDXSRF%NCLN(I)   = Source_NDXSRF%NCLN(I)
   END DO
   DO I = 1,NAN
      Destination_NDXSRF%ATWT(I)    = Source_NDXSRF%ATWT(I)
   END DO
   DO I = 1,NSN
      DO J = 1,4
         Destination_NDXSRF%NDXS(J,I)    = Source_NDXSRF%NDXS(J,I)
      END DO
      DO J = 1,NNS
         Destination_NDXSRF%NOS(J,I)    = Source_NDXSRF%NOS(J,I)
      END DO
      DO J = 1,NON
         Destination_NDXSRF%NOR(J,I)    = Source_NDXSRF%NOR(J,I)
      END DO
   END DO
   DO I = 1,NZONE
      Destination_NDXSRF%VOLZ(I)    = Source_NDXSRF%VOLZ(I)
      Destination_NDXSRF%VFPA(I)    = Source_NDXSRF%VFPA(I)
      Destination_NDXSRF%NSPA(I)    = Source_NDXSRF%NSPA(I)
   END DO
   DO I = 1,NSZ
      Destination_NDXSRF%VLSA(I)    = Source_NDXSRF%VLSA(I)
      Destination_NDXSRF%NSSA(I)    = Source_NDXSRF%NSSA(I)
      Destination_NDXSRF%NZSZ(I)    = Source_NDXSRF%NZSZ(I)
   END DO
END IF

END SUBROUTINE NDXSRF_Copy

!---------------------------------------------------------------------------------------------------------------------------------
!  NDXSRF_DEFINE(USER_NDXSRF,NON,NSN,NNS,NAN,NZONE,NSZ)
!  Allocates and initializes the NDXSRF data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE NDXSRF_DEFINE(USER_NDXSRF,NON,NSN,NNS,NAN,NZONE,NSZ)
IMPLICIT NONE
! Passed in
TYPE (NDXSRF_DATA) USER_NDXSRF  ! A user data variable to be defined by reading in a NDXSRF file
DIF3D_Int :: NON                       ! Number of nuclides in the cross section data
DIF3D_Int :: NSN                       ! Number of nuclide sets (ususually 1)
DIF3D_Int :: NNS                       ! Maximum number of nuclides in a set
DIF3D_Int :: NAN                       ! Number of unique nuclides in the cross section data
DIF3D_Int :: NZONE                     ! Number of zones
DIF3D_Int :: NSZ                       ! Number of subzones
! Local
DIF3D_Int IOS

100 FORMAT('[NDXSRF]...SORRY, BUT I MUST STOP')
101 FORMAT('[NDXSRF]',107('.'))
102 FORMAT('')
105 FORMAT('[NDXSRF]...THERE WAS A FATAL ERROR THAT OCCURED IN (ALLOCATE)',54('.'))

IF (USER_NDXSRF%DEFINED) CALL NDXSRF_VOID(USER_NDXSRF)   ! If the variable is already defined, void it

#ifdef DIF3D_Debug
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF NUCLIDES IN CROSS SECTION DATA",48("."),I16)') NON
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF NUCLIDE SETS IDENTIFIED.......",48("."),I16)') NSN
   WRITE(MODULE_OUT,'("[NDXSRF]...MAXIMUM NUMBER OF NUCLIDES IN ANY SET...",48("."),I16)') NNS
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF DIFFERENT NUCLIDES IN DATA....",48("."),I16)') NAN
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF ZONES.........................",48("."),I16)') NZONE
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF SUBZONES (SUBASSEMBLIES)......",48("."),I16)') NSZ  
#endif

IF (NON .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[NDXSRF]...INVALID NUMBER OF NUCLIDES IN CROSS SECTION DATA",40("."),I16)') NON
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (NSN .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[NDXSRF]...INVALID NUMBER OF NUCLIDE SETS IDENTIFIED.......",40("."),I16)') NSN
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (NNS .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[NDXSRF]...INVALID MAXIMUM NUMBER OF NUCLIDES IN ANY SET...",40("."),I16)') NNS
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (NAN .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[NDXSRF]...INVALID NUMBER OF DIFFERENT NUCLIDES IN DATA....",40("."),I16)') NAN
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (NZONE .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[NDXSRF]...INVALID NUMBER OF ZONES.........................",40("."),I16)') NZONE
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (NSZ .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[NDXSRF]...INVALID NUMBER OF SUBZONES (SUBASSEMBLIES)......",40("."),I16)') NSZ  
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

! Data for card type 0-4
ALLOCATE(USER_NDXSRF%HNNAME(NON),USER_NDXSRF%HANAME(NON),USER_NDXSRF%WPF(NON),USER_NDXSRF%ATWT(NAN),                          &  ! TYPE 2
         USER_NDXSRF%NCLN(NON),USER_NDXSRF%NDXS(4,NSN),USER_NDXSRF%NOS(NNS,NSN),USER_NDXSRF%NOR(NON,NSN),                     &  ! TYPE 2
         USER_NDXSRF%VOLZ(NZONE),USER_NDXSRF%VFPA(NZONE),USER_NDXSRF%VLSA(NSZ),USER_NDXSRF%NSPA(NZONE),USER_NDXSRF%NSSA(NSZ), &  ! TYPE 3
         USER_NDXSRF%NZSZ(NSZ),                                                                                               &  ! TYPE 3
         STAT = IOS)
IF (IOS .NE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[NDXSRF]...FAILED TO ALLOCATE NDXSRF DATA",58("."),I16)') IOS
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF NUCLIDES IN CROSS SECTION DATA",48("."),I16)') NON
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF NUCLIDE SETS IDENTIFIED.......",48("."),I16)') NSN
   WRITE(MODULE_OUT,'("[NDXSRF]...MAXIMUM NUMBER OF NUCLIDES IN ANY SET...",48("."),I16)') NNS
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF DIFFERENT NUCLIDES IN DATA....",48("."),I16)') NAN
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF ZONES.........................",48("."),I16)') NZONE
   WRITE(MODULE_OUT,'("[NDXSRF]...NUMBER OF SUBZONES (SUBASSEMBLIES)......",48("."),I16)') NSZ  
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

USER_NDXSRF%DEFINED  = .TRUE.
USER_NDXSRF%NON   = NON
USER_NDXSRF%NSN   = NSN
USER_NDXSRF%NNS   = NNS
USER_NDXSRF%NAN   = NAN
USER_NDXSRF%NZONE = NZONE
USER_NDXSRF%NSZ   = NSZ
! Initialize arrays
USER_NDXSRF%HNNAME = 'UNDEFINE'
USER_NDXSRF%HANAME = 'UNDEFINE'
USER_NDXSRF%HNNAME = 'UNDEFINE'
USER_NDXSRF%WPF    = 0.0
USER_NDXSRF%ATWT   = 0.0
USER_NDXSRF%NCLN   = 0
USER_NDXSRF%NDXS   = 0
USER_NDXSRF%NOS    = 0
USER_NDXSRF%NOR    = 0
USER_NDXSRF%VOLZ   = 0.0
USER_NDXSRF%VFPA   = 0.0
IF (NSZ .GT. 0) USER_NDXSRF%VLSA   = 0.0
USER_NDXSRF%NSPA   = 0
IF (NSZ .GT. 0) USER_NDXSRF%NSSA   = 0
IF (NSZ .GT. 0) USER_NDXSRF%NZSZ   = 0

END SUBROUTINE NDXSRF_DEFINE

!---------------------------------------------------------------------------------------------------------------------------------
!  NDXSRF_VOID(USER_NDXSRF)
!  Deallocates the data
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE NDXSRF_VOID(USER_NDXSRF)
IMPLICIT NONE
TYPE (NDXSRF_DATA) USER_NDXSRF  ! A user data variable to be defined by reading in a NDXSRF file
! local
DIF3D_Int IOS

100 FORMAT('[NDXSRF]...SORRY, BUT I MUST STOP')
101 FORMAT('[NDXSRF]',107('.'))
102 FORMAT('')
105 FORMAT('[NDXSRF]...THERE WAS A FATAL ERROR THAT OCCURED IN (VOIDTYPE)',54('.'))

IF (USER_NDXSRF%DEFINED) THEN
   DEALLOCATE(USER_NDXSRF%HNNAME,USER_NDXSRF%HANAME,USER_NDXSRF%WPF,USER_NDXSRF%ATWT,                    &  ! TYPE 2
              USER_NDXSRF%NCLN,USER_NDXSRF%NDXS,USER_NDXSRF%NOS,USER_NDXSRF%NOR,                         &  ! TYPE 2
              USER_NDXSRF%VOLZ,USER_NDXSRF%VFPA,USER_NDXSRF%VLSA,USER_NDXSRF%NSPA,USER_NDXSRF%NSSA,      &  ! TYPE 3
              USER_NDXSRF%NZSZ,                                                                          &  ! TYPE 3
              STAT = IOS)
   IF (IOS .NE. 0) THEN
      WRITE(MODULE_OUT,105)
      WRITE(MODULE_OUT,'("[NDXSRF]...FAILED TO DEALLOCATE NDXSRF DATA",72("."))')
      WRITE(MODULE_OUT,100)
      CALL Abort
   END IF
   USER_NDXSRF%DEFINED = .FALSE.
END IF

END SUBROUTINE NDXSRF_VOID

END MODULE NDXSRF_io
