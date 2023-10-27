!---------------------------------------------------------------------------------------------------------------------------------
!  Using the historical format for the COMPXS file, import the data from the binary file 
!  This version assumes the entire data set can fit into memory with full block scattering
!---------------------------------------------------------------------------------------------------------------------------------
!  COMPXS_ASSIGNPRINTINFO(OUTPUT_UNIT)                                  !  Sets the output unit and DEBUG level settings for the module
!  COMPXS_VOID(USER_COMPXS)                                             !  Provides a path to void the COMPXS data type
!  COMPXS_DEFINE(USER_COMPXS,NCMP,NGROUP,ISCHI,NFAM,MAXORD)             !  Defines a COMPXS data structure based upon user input
!  COMPXS_COPY(SOURCE_COMPXS,DESTINATION_COMPXS)                        !  Copies the data from SOURCE_COMPXS to DESTINATION_COMPXS
!  COMPXS_ADDSCALE(TargetFactor,Target_COMPXS,AddOnFactor,AddOn_COMPXS) !  Scales the data in Target_COMPXS and adds on AddOn_COMPXS with a factor
!  COMPXS_ZERO(USER_COMPXS)                                             !  Zeros all of the arrays in the COMPXS data structure
!  COMPXS_PRINT(USER_COMPXS)                                            !  Prints a COMPXS data structure
!  COMPXS_BALANCEXS(USER_COMPXS,L_USETOTAL)                             !  Returns a DIF3D_Log variable describing the success, if any (L_USETOTAL) is used to switch TOTAL from default of TRANSPORT
!---------------------------------------------------------------------------------------------------------------------------------
!  COMPXS_UPDATEPROPERTIES(USER_COMPXS)                                 !  Update several property variable pieces of the existing COMPXS data set before exporting it
!---------------------------------------------------------------------------------------------------------------------------------
MODULE COMPXS_io 
#include "DIF3D_Types.h"
IMPLICIT NONE

! This data type contains all of the scattering data in full block form
TYPE COMPXS_DATA
   DIF3D_Log :: DEFINED=.FALSE.
   DIF3D_Log :: ARRAYSDEFINED = .FALSE.
   DIF3D_FileNames FILENAME                                          ! Originating file name
   ! CARD TYPE 1
   ! NCMP,NGROUP,ISCHI,NFCMP,MAXUP,MAXDN,NFAM,MAXORD,NDUM2,NDUM3
   COMPXS_Int :: NCMP                                     ! NUMBER OF COMPOSITIONS.
   COMPXS_Int :: NGROUP                                   ! NUMBER OF ENERGY GROUPS.
   COMPXS_Int :: ISCHI                                    ! PROMPT FISSION SPECTRUM FLAG. 
                                                                  !   ISCHI=0      IF THERE IS NO SET-WIDE PROMPT CHI. 
                                                                  !   ISCHI=1      IF THERE IS A SET-WIDE PROMPT CHI VECTOR.
                                                                  !   ISCHI=NGROUP IF THERE IS A SET-WIDE PROMPT CHI MATRIX.
   COMPXS_Int :: NFCMP                                    ! NUMBER OF FISSIONABLE COMPOSITIONS.
   COMPXS_Int :: MAXUP                                    ! MAXIMUM NUMBER OF GROUPS OF UPSCATTERING FOR THE SET.
   COMPXS_Int :: MAXDN                                    ! MAXIMUM NUMBER OF GROUPS OF DOWNSCATTERING FOR THE SET.
   COMPXS_Int :: NFAM                                     ! NUMBER OF DELAYED NEUTRON FAMILIES.
   COMPXS_Int :: MAXORD                                   ! ANISOTROPIC SCATTERING ORDER
   COMPXS_Int :: LEGTERMS                                 ! NUMBER OF LEGENDRE TERMS ! MAXORD + 1
   COMPXS_Int :: NDUM2                                    ! RESERVED.
   COMPXS_Int :: NDUM3                                    ! RESERVED.
   ! CARD TYPE 2   COMPOSITION INDEPENDENT DATA
   COMPXS_Real,DIMENSION(:,:),    POINTER :: CHI       ! ISCHI,NGROUP ; PROMPT FISSION FRACTION INTO GROUP J FROM GROUP I. 
   !                                                                             ; IF ISCHI=1, THE LIST REDUCES TO (CHI(J),J=1,NGROUP), 
   !                                                                             ; WHERE CHI(J) IS THE FISSION FRACTION INTO GROUP J.
   COMPXS_Real,DIMENSION(:),      POINTER :: VEL       ! NGROUP       ; MEAN NEUTRON VELOCITY IN GROUP J (CM/SEC).
   COMPXS_Real,DIMENSION(:),      POINTER :: EMAX      ! NGROUP+1     ; MAXIMUM ENERGY BOUND OF GROUP J (EV).
   COMPXS_Real,DIMENSION(:,:),    POINTER :: CHID      ! NGROUP,NFAM  ; FRACTION OF DELAYED NEUTRONS EMITTED INTO NEUTRON ENERGY GROUP J FROM PRECURSOR FAMILY K.
   COMPXS_Real,DIMENSION(:),      POINTER :: FLAM      ! NFAM         ; DELAYED NEUTRON PRECURSOR DECAY CONSTANT FOR FAMILY K.
   COMPXS_Int,DIMENSION(:),      POINTER :: NKFAM     ! NCMP         ; NUMBER OF FAMILIES TO WHICH FISSION IN COMPOSITION J CONTRIBUTES DELAYED NEUTRON PRECURSORS.
   ! CARD TYPE 3   COMPOSITION SPECIFICATIONS                  
   COMPXS_Int,DIMENSION(:),      POINTER :: ICHI      ! NCMP         ; PROMPT FISSION SPECTRUM FLAG FOR THIS COMPOSITION. 
   !                                                                             ; ICHI=-1      IF COMPOSITION USES THE SET-WIDE PROMPT CHI GIVEN IN SET CHI RECORD (BELOW). 
   !                                                                             ; ICHI=0       IF COMPOSITION IS NOT FISSIONABLE. 
   !                                                                             ; ICHI=1       FOR COMPOSITION PROMPT CHI VECTOR. 
   !                                                                             ; ICHI=NGROUP FOR COMPOSITION PROMPT CHI MATRIX.
   COMPXS_Int,DIMENSION(:,:),    POINTER :: NUP       ! NGROUP,NCMP  ; NUMBER OF GROUPS OF UPSCATTERING INTO GROUP I FROM LOWER ENERGY GROUPS FOR THE CURRENT COMPOSITION
   COMPXS_Int,DIMENSION(:,:),    POINTER :: NDN       ! NGROUP,NCMP  ; NUMBER OF GROUPS OF DOWNSCATTERING INTO GROUP I FROM HIGHER ENERGY GROUPS FOR THE CURRENT COMPOSITION
   COMPXS_Int,DIMENSION(:,:),    POINTER :: NUMFAM    ! NFAM,NCMP    ; FAMILY NUMBER OF THE I-TH YIELD VECTOR IN ARRAY SNUDEL(I).
   ! CARD TYPE 4   COMPOSITION MACROSCO PIC GROUP CROSS SECTIONS
   COMPXS_Real,DIMENSION(:,:),    POINTER :: XA        ! NGROUP,NCMP                 ; ABSORPTION CROSS SECTION.
   COMPXS_Real,DIMENSION(:,:),    POINTER :: XTOT      ! NGROUP,NCMP                 ; TOTAL CROSS SECTION.
   COMPXS_Real,DIMENSION(:,:),    POINTER :: XREM      ! NGROUP,NCMP                 ; REMOVAL CROSS SECTION, TOTAL CROSS SECTION FOR REMOVING A NEUTRON FROM GROUP J DUE TO ALL PROCESSES.
   COMPXS_Real,DIMENSION(:,:),    POINTER :: XTR       ! NGROUP,NCMP                 ; TRANSPORT CROSS SECTION.
   COMPXS_Real,DIMENSION(:,:),    POINTER :: XF        ! NGROUP,NCMP                 ; FISSION CROSS SECTION, PRESENT ONLY IF ICHI.NE.0. 
   COMPXS_Real,DIMENSION(:,:),    POINTER :: XNF       ! NGROUP,NCMP                 ; TOTAL NUMBER OF NEUTRONS EMITTED PER FISSION TIMES XF, PRESENT ONLY IF ICHI.NE.0.
   COMPXS_Real,DIMENSION(:,:,:),  POINTER :: ZONE_CHI  ! NGROUP,NGROUP,NCMP          ; PROMPT FISSION FRACTION INTO GROUP J FROM GROUP I, PRESENT ONLY IF ICHI.GT.0. 
   !                                                     J      ICHI                 ;    IF ICHI=1,THE LIST REDUCES TO THE SINGLE NUMBER CHI, WHICH IS THE PROMPT FISSION FRACTION INTO GROUP J.
   COMPXS_Real,DIMENSION(:,:,:,:),POINTER :: XSCAT     ! NGROUP,NGROUP,LEGTERMS,NCMP ; TOTAL SCATTERING CROSS SECTION
                                                       ! Stored as scatter from (Group->,Group) not (Group,<-Group), sorry...MAS
   COMPXS_Real,DIMENSION(:,:),    POINTER :: PC        ! NGROUP,NCMP                 ; PC TIMES THE GROUP J REGION INTEGRATED FLUX FOR THE 
   !                                                                                            ;  REGIONS CONTAINING THE CURRENT COMPOSITION YIELDS 
   !                                                                                            ;  THE POWER IN WATTS IN THOSE REGIONS AND ENERGY GROUP
   !                                                                                            ;  J DUE TO FISSIONS AND NON-FISSION ABSORPTIONS.
   COMPXS_Real,DIMENSION(:,:),    POINTER :: A1        ! NGROUP,NCMP                 ; FIRST DIMENSION DIRECTIONAL DIFFUSION COEFFICIENT MULTIPLIER.
   COMPXS_Real,DIMENSION(:,:),    POINTER :: B1        ! NGROUP,NCMP                 ; FIRST DIMENSION DIRECTIONAL DIFFUSION COEFFICIENT ADDITIVE TERM.
   COMPXS_Real,DIMENSION(:,:),    POINTER :: A2        ! NGROUP,NCMP                 ; SECOND DIMENSION DIRECTIONAL DIFFUSION COEFFICIENT MULTIPLIER.
   COMPXS_Real,DIMENSION(:,:),    POINTER :: B2        ! NGROUP,NCMP                 ; SECOND DIMENSION DIRECTIONAL DIFFUSION COEFFICIENT ADDITIVE TERM.
   COMPXS_Real,DIMENSION(:,:),    POINTER :: A3        ! NGROUP,NCMP                 ; THIRD DIMENSION DIRECTIONAL DIFFUSION COEFFICIENT MULTIPLIER.
   COMPXS_Real,DIMENSION(:,:),    POINTER :: B3        ! NGROUP,NCMP                 ; THIRD DIMENSION DIRECTIONAL DIFFUSION COEFFICIENT ADDITIVE TERM.
   COMPXS_Real,DIMENSION(:,:,:),  POINTER :: SNUDEL    ! NFAM,NGROUP,NCMP            ; NUMBER OF DELAYED NEUTRON PRECURSORS PRODUCED IN FAMILY NUMBER NUMFAM(I) PER FISSION IN GROUP J.
   COMPXS_Real,DIMENSION(:,:),    POINTER :: XN2N      ! NGROUP,NCMP                 ; N,2N REACTION CROSS SECTION
   !                                                              ! THE MACROSCOPIC XN2N(J) TIMES THE FLUX IN GROUP J GIVES THE RATE AT WHICH N,2N REACTIONS OCCUR IN GROUP J.  
   !                                                              ! THUS, FOR N,2N SCATTERING, XN2N(J)=0.5*(SUM OF SCAT(J TO G)) SUMMED OVER ALL G WHERE SCAT IS THE N,2N SCATTERING MATRIX.
   ! CARD TYPE 5   POWER CONVERSION FACTORS
   COMPXS_Real,DIMENSION(:),      POINTER :: FPWS      ! NCMP                        ; FISSIONS/WATT-SECOND FOR EACH COMPOSITION
   COMPXS_Real,DIMENSION(:),      POINTER :: CPWS      ! NCMP                        ; CAPTURES/WATT-SECOND FOR EACH COMPOSITION
END TYPE

TYPE (COMPXS_DATA), SAVE :: MODULE_COMPXS

DIF3D_Int :: MODULE_OUT = 6                                                          ! Module output unit
PRIVATE MODULE_OUT                                                                   ! Variables

CONTAINS

!---------------------------------------------------------------------------------------------------------------------------------
!  COMPXS_ASSIGNPRINTINFO(OUTPUT_UNIT,DEBUGPRINT)
!  Sets the output unit and DEBUG level settings for the module
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE COMPXS_ASSIGNPRINTINFO(MODULE_OUTPUT_UNIT)
DIF3D_Int MODULE_OUTPUT_UNIT    ! Reset the object output unit

MODULE_OUT = MODULE_OUTPUT_UNIT

END SUBROUTINE COMPXS_ASSIGNPRINTINFO

!---------------------------------------------------------------------------------------------------------------------------------
!  COMPXS_PRINT(USER_COMPXS)
!  Prints a COMPXS data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE COMPXS_PRINT(USER_COMPXS)
IMPLICIT NONE
! Passed in
TYPE (COMPXS_DATA) USER_COMPXS                   ! A user data variable to be defined by reading in a COMPXS file
! Local
DIF3D_Int IOS,G,I,J,K,L,M,N,LLL
DIF3D_Int CURC,SCATLENGTH,SCATSTART
DIF3D_Int INPUTUNIT
DIF3D_Int, PARAMETER :: LocalMaxRemoval=8
DIF3D_Real SCATREMOVAL(LocalMaxRemoval) ! MAS added this to print out the even- and odd- parity removal cross sections.

100 FORMAT('[COMPXS]...SORRY, BUT I MUST STOP')
101 FORMAT('[COMPXS]',107('.'))
102 FORMAT('')
105 FORMAT('[COMPXS]...THERE WAS A FATAL ERROR THAT OCCURED IN (PRINT)',57('.'))

IF (.NOT. USER_COMPXS%DEFINED) THEN
   WRITE(MODULE_OUT,101)
   WRITE(MODULE_OUT,'("[COMPXS]...USER COMPXS VARIABLE NOT DEFINED FOR PRINTING...NONFATAL",41("."))')
   WRITE(MODULE_OUT,101)
ELSE
   WRITE(MODULE_OUT,101)
   WRITE(MODULE_OUT,'("[COMPXS]...NUMBER OF COMPOSITIONS............",54("."),I16)') USER_COMPXS%NCMP
   WRITE(MODULE_OUT,'("[COMPXS]...NUMBER OF GROUPS..................",54("."),I16)') USER_COMPXS%NGROUP
   WRITE(MODULE_OUT,'("[COMPXS]...PROMPT FISSION SPECTRUM FLAG......",54("."),I16)') USER_COMPXS%ISCHI
   WRITE(MODULE_OUT,'("[COMPXS]...NUMBER OF FISSIONABLE COMPOSITIONS",54("."),I16)') USER_COMPXS%NFCMP
   WRITE(MODULE_OUT,'("[COMPXS]...MAXIMUM UPSCATTERING GROUPS.......",54("."),I16)') USER_COMPXS%MAXUP
   WRITE(MODULE_OUT,'("[COMPXS]...MAXIMUM DOWNSCATTERING GROUPS.....",54("."),I16)') USER_COMPXS%MAXDN
   WRITE(MODULE_OUT,'("[COMPXS]...NUMBER OF DELAYED NEUTRON FAMILIES",54("."),I16)') USER_COMPXS%NFAM
   WRITE(MODULE_OUT,'("[COMPXS]...ANISOTROPIC SCATTERING ORDER......",54("."),I16)') USER_COMPXS%MAXORD
   WRITE(MODULE_OUT,'("[COMPXS]...NUMBER OF Legendre Terms..........",54("."),I16)') USER_COMPXS%LEGTERMS
   WRITE(MODULE_OUT,'("[COMPXS]...RESERVED 1?.......................",54("."),I16)') USER_COMPXS%NDUM2
   WRITE(MODULE_OUT,'("[COMPXS]...RESERVED 2?.......................",54("."),I16)') USER_COMPXS%NDUM3

   LLL = USER_COMPXS%LEGTERMS
   IF (LLL .GT. LocalMaxRemoval) LLL = LocalMaxRemoval

   WRITE(MODULE_OUT,'("[COMPXS]... GROUP    LOWER E BOUND    UPPER E BOUND      VELOCITY   ",47("."))')
   DO G = 1,USER_COMPXS%NGROUP
      WRITE(MODULE_OUT,'("[COMPXS]...",I6,1X,1PE16.7,1X,1PE16.7,1X,1PE16.7,47("."))') &
                  G,USER_COMPXS%EMAX(G+1),USER_COMPXS%EMAX(G),USER_COMPXS%VEL(G)
   END DO
   IF (USER_COMPXS%ISCHI .GT. 0) THEN
      WRITE(MODULE_OUT,'("[COMPXS]...FILE WIDE CHI DATA",86("."))')
      WRITE(MODULE_OUT,'("[COMPXS]...",6X,1000(5X,"CHI ",I4))') (I,I=1,USER_COMPXS%ISCHI)
      DO G = 1,USER_COMPXS%NGROUP
         WRITE(MODULE_OUT,'("[COMPXS]...",I4,2X,1000(1PE13.6))') G,(USER_COMPXS%CHI(I,G),I=1,USER_COMPXS%ISCHI)
      END DO
   END IF
IF (USER_COMPXS%NFAM .GT. 0) THEN
   DO G = 1,USER_COMPXS%NGROUP
      WRITE(MODULE_OUT,'("[COMPXS]...FRACTION OF DELAYED NEUTRONS FROM EACH PRECURSOR FAMILY EMITTED INTO GROUP",14("."),I16)') G
      WRITE(MODULE_OUT,'("[COMPXS]...",6(1X,I3,1P1E12.3))')(K,USER_COMPXS%CHID(G,K),K=1,USER_COMPXS%NFAM)
   END DO
   WRITE(MODULE_OUT,'("[COMPXS]...DELAYED NEUTRON PRECURSOR DECAY CONSTANT FOR EACH FAMILY",48("."))')
   WRITE(MODULE_OUT,'("[COMPXS]...",6(1X,I3,1P1E13.6))')(K,USER_COMPXS%FLAM(K),K=1,USER_COMPXS%NFAM)
   WRITE(MODULE_OUT,'("[COMPXS]...NUMBER OF FAMILIES TO WHICH COMPOSITION FISSION CONTRIB. DELAYED NEUTRON PRECURSORS",21("."))')
   WRITE(MODULE_OUT,'("[COMPXS]...",6(1X,I3,I12))')(K,USER_COMPXS%NKFAM(K),K=1,USER_COMPXS%NCMP)
END IF

DO CURC = 1,USER_COMPXS%NCMP
      WRITE (MODULE_OUT,'("[COMPXS]...CROSS SECTION DATA FOR COMPOSITION ORDERED",46("."),I16)') CURC
      WRITE (MODULE_OUT,1100)
                     1100 FORMAT("[COMPXS]...",3X,'GROUP',  3X,'TRANS (P1)',3X,'TOTAL (P0)',4X,'ABSORP.',      &
                                 6X,'REMOVAL',8X,'N2N', 8X,'FISSION',4X,'YIELD/FISS',2x,'....')
      DO G = 1,USER_COMPXS%NGROUP
         WRITE(MODULE_OUT,'("[COMPXS]...",2X,I6,7(1PE13.5),1X,4("."))')G,USER_COMPXS%XTR(G,CURC),USER_COMPXS%XTOT(G,CURC),  &
                                               USER_COMPXS%XA(G,CURC),USER_COMPXS%XREM(G,CURC),USER_COMPXS%XN2N(G,CURC),   &
                                                USER_COMPXS%XF(G,CURC),USER_COMPXS%XNF(G,CURC)
      END DO
      WRITE(MODULE_OUT,'("[COMPXS]...ICHI = ",I6," FOR COMPOSITION",59("."),I16)') USER_COMPXS%ICHI(CURC),CURC
      IF (USER_COMPXS%ICHI(CURC) .GT. 0) THEN
         WRITE (MODULE_OUT,1110) (I,I=1,USER_COMPXS%ICHI(CURC))
                            1110 FORMAT("[COMPXS]...",8X,1000(5X,'CHI ',I3)) 
         DO G = 1,USER_COMPXS%NGROUP
            WRITE(MODULE_OUT,'("[COMPXS]...",2X,I6,1P,1000E12.5)') G,(USER_COMPXS%ZONE_CHI(G,I,CURC),I=1,USER_COMPXS%ICHI(CURC))
         END DO
      END IF

      WRITE(MODULE_OUT,1120)
                    1120 FORMAT("[COMPXS]...",3X,'GROUP',  2X,'POW CONVER.',  4X,'DIFF A1',  6X,'DIFF B1',     &
                                6X,'DIFF A2',6X,'DIFF B2',6X,'DIFF A3',6X,'DIFF B3',3X,4("."))
      DO G=1,USER_COMPXS%NGROUP
         WRITE(MODULE_OUT,'("[COMPXS]...",2X,I6,7(1PE13.5),1x,4("."))')G,USER_COMPXS%PC(G,CURC),USER_COMPXS%A1(G,CURC),  &
                   USER_COMPXS%B1(G,CURC),USER_COMPXS%A2(G,CURC),USER_COMPXS%B2(G,CURC),USER_COMPXS%A3(G,CURC),            &
                   USER_COMPXS%B3(G,CURC)
      END DO

      WRITE(MODULE_OUT,1121) (G-1,G=1,8)
      1121 FORMAT('[COMPXS]...',' GROUP ',10(2X,'P',I1,' Removal',1X))
      DO G=1,USER_COMPXS%NGROUP
         DO L = 1,LocalMaxRemoval
            SCATREMOVAL(L) = USER_COMPXS%XTOT(G,CURC)
         END DO
         DO L = 1,LLL
            SCATREMOVAL(L) = SCATREMOVAL(L) - USER_COMPXS%XSCAT(G,G,L,CURC)
         END DO
         WRITE(MODULE_OUT,'("[COMPXS]...",1X,I5,1X,10(1PE13.5))') G,(SCATREMOVAL(L),L=1,LocalMaxRemoval)
      END DO

      DO I = 1,USER_COMPXS%NKFAM(CURC)
         WRITE(MODULE_OUT,                                                                                                  &
 '("[COMPXS]...NUMBER OF DELAYED NEUTRON PRECURSORS PRODUCED IN FAMILY NUMBER ",I3," PER FISSION IN EACH GROUP",12("."))')  &
               USER_COMPXS%NUMFAM(I,CURC)
         WRITE(MODULE_OUT,'("[COMPXS]...",6(1X,I3,1PE13.6))') (G,USER_COMPXS%SNUDEL(I,G,CURC),G=1,USER_COMPXS%NGROUP)
      END DO

      DO L = 1,USER_COMPXS%LEGTERMS
         WRITE(MODULE_OUT,'("[COMPXS]...SCATTERING MATRIX FOR LEGENDRE ORDER",52("."),I16)') L-1
         WRITE(MODULE_OUT,'("[COMPXS]...",4X,1000(2X,"TO GROUP ",I3))') (I,I=1,USER_COMPXS%NGROUP)  ! Yes dumbass, it is (Group->,Group) not (Group,<-Group), Sorry...MAS
         DO G=1,USER_COMPXS%NGROUP
            WRITE(MODULE_OUT,'("[COMPXS]...",I4,1000(1PE13.5,1X))') G,(USER_COMPXS%XSCAT(G,I,L,CURC),I=1,USER_COMPXS%NGROUP)
         END DO
      END DO
END DO ! CURRENTCOMP
! CARD TYPE 5   POWER CONVERSION FACTORS
   WRITE(MODULE_OUT,'("[COMPXS]...FISSION POWER CONVERSION FACTOR",73("."))')
   WRITE(MODULE_OUT,'(("[COMPXS]...",6(1X,I3,1PE13.6)))')(I,USER_COMPXS%FPWS(I),I=1,USER_COMPXS%NCMP)
   WRITE(MODULE_OUT,'("[COMPXS]...CAPTURE POWER CONVERSION FACTOR",73("."))')
   WRITE(MODULE_OUT,'(("[COMPXS]...",6(1X,I3,1PE13.6)))')(I,USER_COMPXS%CPWS(I),I=1,USER_COMPXS%NCMP)
END IF

END SUBROUTINE COMPXS_PRINT

!---------------------------------------------------------------------------------------------------------------------------------
!  COMPXS_BALANCEXS(USER_COMPXS,L_USETOTAL)
!  Returns a DIF3D_Log variable describing the success, if any (L_USETOTAL) is used to switch TOTAL from default of TRANSPORT
! 
!  Balance the XS using the following method where Total and absorption are assumed to be accurate.
!  Total = Absorption + Scattering (Adjust within group scattering)
!  Removal = Absorption + Scattering - within group scattering = Total - new within group scattering
!  Transport = Transport     no hope to fix here.
!  Fission = Fission but clear out small negative values and issue warning
!  NuFission = NuFission but clear out small negative values and issue warning
!  Chi       = Normalized Chi
!  
!---------------------------------------------------------------------------------------------------------------------------------
DIF3D_Log FUNCTION COMPXS_BALANCEXS(USER_COMPXS,L_USETOTAL)
IMPLICIT NONE
! Passed in
TYPE (COMPXS_DATA) USER_COMPXS                   ! A user data variable to be defined by reading in a COMPXS file
DIF3D_Log L_USETOTAL
! Local variables
DIF3D_Log L_BALANCEERROR
DIF3D_Int G,GG,CURC
REAL CONST1,CONST2

100 FORMAT('[COMPXS]...SORRY, BUT I MUST STOP')
101 FORMAT('[COMPXS]',107('.'))
102 FORMAT('')
105 FORMAT('[COMPXS]...THERE WAS A BALANCE ERROR THAT OCCURED IN (BALANCEXS)',51('.'))

L_BALANCEERROR = .FALSE.
! Sweep through all compositions and balance cross sections
DO CURC = 1,USER_COMPXS%NCMP
   DO G = 1,USER_COMPXS%NGROUP
      ! Adjust within group scattering such that transport cross section is total (i.e. setup for a diffusion calculation)
      CONST1 = 0.0                              ! out of group scattering portion
      DO GG = 1,USER_COMPXS%NGROUP
         IF (G .NE. GG) CONST1 = CONST1 + USER_COMPXS%XSCAT(G,GG,1,CURC)
      END DO
      IF (L_USETOTAL) THEN
         CONST2 = USER_COMPXS%XTOT(G,CURC) - CONST1 - USER_COMPXS%XA(G,CURC)
      ELSE
         CONST2 = USER_COMPXS%XTR(G,CURC) - CONST1 - USER_COMPXS%XA(G,CURC)
      END IF
      IF (CONST2 .LT. 0.0D0) THEN
         WRITE(MODULE_OUT,105)
         WRITE(MODULE_OUT,200) CURC,G,CONST2
         200 FORMAT('[COMPXS]...NEGATIVE ISOTROPIC SCATTERING FOR COMPOSITION ',I5,' GROUP ',I6,24("."),1PE16.8)
         L_BALANCEERROR = .TRUE.
      ELSE
         USER_COMPXS%XSCAT(G,G,1,CURC) = CONST2
      END IF
      ! Adjust removal
      IF (L_USETOTAL) THEN
         CONST2 = USER_COMPXS%XTOT(G,CURC) - USER_COMPXS%XSCAT(G,G,1,CURC)
      ELSE
         CONST2 = USER_COMPXS%XTR(G,CURC) - USER_COMPXS%XSCAT(G,G,1,CURC)
      END IF

      IF (CONST2 .LT. 0.0D0) THEN
         WRITE(MODULE_OUT,105)
         WRITE(MODULE_OUT,210) CURC,G,CONST2
         210 FORMAT('[COMPXS]...NEGATIVE REMOVAL FOR COMPOSITION ',I5,' GROUP ',I6,37("."),1PE16.8)
         L_BALANCEERROR = .TRUE.
      ELSE
         USER_COMPXS%XREM(G,CURC) = CONST2
      END IF
      ! Check transport XS (in general this should be similar to total)
      ! Transport in general = 1 / (Scattering * ( 1- (2/3A)))
      IF ((USER_COMPXS%XTR(G,CURC) .LT. 0.75d0*USER_COMPXS%XTOT(G,CURC)) .OR.  &
          (USER_COMPXS%XTR(G,CURC) .GT. 1.25d0*USER_COMPXS%XTOT(G,CURC))        ) THEN
         WRITE(MODULE_OUT,220) CURC,G,USER_COMPXS%XTR(G,CURC)
         220 FORMAT('[COMPXS]...SIGNIFICANT TRANS/TOTAL DISAGREEMENT FOR COMP. ',I5,' GROUP ',I6,   &
                     ' TRANS=',1PE12.5,' TOTAL= ',1PE12.5,'.'                                                )
      END IF
      ! Check and adjust fission XS
      IF (USER_COMPXS%XF(G,CURC) .LT. 0.0D0) THEN
         WRITE(MODULE_OUT,230) CURC,G,USER_COMPXS%XF(G,CURC)
         230 FORMAT('[COMPXS]...FISSION XS IS NEGATIVE FOR COMPOSITION ',I5,' GROUP ',I6,31("."),1PE16.8)
         L_BALANCEERROR = .TRUE.
         USER_COMPXS%XF(G,CURC) = 0.0D0
      END IF
      ! Check and adjust nu-fission XS
      IF (USER_COMPXS%XNF(G,CURC) .LT. 0.0D0) THEN
         WRITE(MODULE_OUT,240) CURC,G,USER_COMPXS%XNF(G,CURC)
         240 FORMAT('[COMPXS]...NU-FISSION XS IS NEGATIVE FOR COMPOSITION ',I5,' GROUP ',I6,28("."),1PE16.8)
         L_BALANCEERROR = .TRUE.
         USER_COMPXS%XNF(G,CURC) = 0.0D0
      END IF
      ! Check and adjust Power conversion
      IF (USER_COMPXS%PC(G,CURC) .LT. 0.0D0) THEN
         WRITE(MODULE_OUT,250) CURC,G,USER_COMPXS%PC(G,CURC)
         250 FORMAT('[COMPXS]...POWER CONVERSION IS NEGATIVE FOR COMPOSITION ',I5,' GROUP ',I6,25("."),1PE16.8)
         L_BALANCEERROR = .TRUE.
         USER_COMPXS%PC(G,CURC) = 0.0D0
      END IF
   END DO
END DO

COMPXS_BALANCEXS = L_BALANCEERROR

END FUNCTION COMPXS_BALANCEXS

!---------------------------------------------------------------------------------------------------------------------------------
!  COMPXS_DEFINE(USER_COMPXS,NCMP,NGROUP,ISCHI,NFAM,MAXORD)
!  Allocates and initializes the data structures in this module based on user info
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE COMPXS_DEFINE(USER_COMPXS,NCMP,NGROUP,ISCHI,NFAM,MAXORD)
IMPLICIT NONE
! Passed in
TYPE (COMPXS_DATA) USER_COMPXS                   ! A user data variable to be defined by reading in a COMPXS file
DIF3D_Int NCMP         ! Number of compositions
DIF3D_Int NGROUP       ! Number of groups
DIF3D_Int ISCHI        ! Chi type, 1 or NGROUP
DIF3D_Int NFAM         ! Number of families
DIF3D_Int MAXORD       ! Number of Legendre expansion terms
! Local
DIF3D_Int IOS,LEGTERMS

100 FORMAT('[COMPXS]...SORRY, BUT I MUST STOP')
101 FORMAT('[COMPXS]',107('.'))
102 FORMAT('')
105 FORMAT('[COMPXS]...THERE WAS A FATAL ERROR THAT OCCURED IN (ALLOCATE)',54('.'))

IF ((USER_COMPXS%DEFINED) .OR. (USER_COMPXS%ARRAYSDEFINED)) CALL COMPXS_VOID(USER_COMPXS)

! Check the variables to ensure sanity
IF (NCMP .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[COMPXS]...INVALID NUMBER OF COMPOSITIONS",58("."),I16)') NCMP
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (NGROUP .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[COMPXS]...INVALID NUMBER OF GROUPS",64("."),I16)') NGROUP
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (NFAM .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[COMPXS]...INVALID NUMBER OF FAMILIES",62("."),I16)') NFAM
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

IF (MAXORD .LT. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[COMPXS]...INVALID ANISOTROPIC SCATTERING ORDER",52("."),I16)') MAXORD
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

LEGTERMS = MAXORD + 1

IF ((ISCHI .NE. 0) .AND. (ISCHI .NE. 1) .AND. (ISCHI .NE. NGROUP)) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[COMPXS]...INVALID NUMBER OF CHI TERMS (1 OR NGROUP)",47("."),I16)') ISCHI
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

#ifdef DIF3D_Debug
   WRITE(MODULE_OUT,'("[COMPXS]...NUMBER OF COMPOSITIONS",      66("."),I16)') NCMP
   WRITE(MODULE_OUT,'("[COMPXS]...NUMBER OF ENERGY GROUPS",     65("."),I16)') NGROUP
   WRITE(MODULE_OUT,'("[COMPXS]...PROMPT FISSION SPECTRUM FLAG",60("."),I16)') ISCHI
   WRITE(MODULE_OUT,'("[COMPXS]...ANISOTROPIC SCATTERING ORDER",60("."),I16)') MAXORD
   WRITE(MODULE_OUT,'("[COMPXS]...NUMBER OF FAMILIES",          70("."),I16)') NFAM
#endif

! Allocate the data arrays
ALLOCATE(USER_COMPXS%VEL(NGROUP),USER_COMPXS%EMAX(NGROUP+1),                                                                     &
         USER_COMPXS%CHID(NGROUP,NFAM),USER_COMPXS%FLAM(NFAM),USER_COMPXS%NKFAM(NCMP),                                           &
         USER_COMPXS%NUMFAM(NFAM,NCMP),USER_COMPXS%SNUDEL(NFAM,NGROUP,NCMP),                                                     &
         USER_COMPXS%CHI(ISCHI,NGROUP), USER_COMPXS%ICHI(NCMP),USER_COMPXS%ZONE_CHI(NGROUP,NGROUP,NCMP),                         &
         USER_COMPXS%NUP(NGROUP,NCMP), USER_COMPXS%NDN(NGROUP,NCMP),                                                             &
         USER_COMPXS%XA(NGROUP,NCMP),USER_COMPXS%XTOT(NGROUP,NCMP),USER_COMPXS%XREM(NGROUP,NCMP),USER_COMPXS%XTR(NGROUP,NCMP),   &
         USER_COMPXS%XF(NGROUP,NCMP),USER_COMPXS%XNF(NGROUP,NCMP),USER_COMPXS%XSCAT(NGROUP,NGROUP,LEGTERMS,NCMP),                &
         USER_COMPXS%XN2N(NGROUP,NCMP),USER_COMPXS%PC(NGROUP,NCMP),                                                              &
         USER_COMPXS%A1(NGROUP,NCMP),USER_COMPXS%A2(NGROUP,NCMP),USER_COMPXS%A3(NGROUP,NCMP),                                    &
         USER_COMPXS%B1(NGROUP,NCMP),USER_COMPXS%B2(NGROUP,NCMP),USER_COMPXS%B3(NGROUP,NCMP),                                    &
         USER_COMPXS%FPWS(NCMP),USER_COMPXS%CPWS(NCMP),                                                                          &
         STAT = IOS                                                                                                              )
IF (IOS .NE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[COMPXS]...FAILED TO ALLOCATE COMPXS DATA",74("."))')
   WRITE(MODULE_OUT,100)
   CALL Abort
END IF

USER_COMPXS%DEFINED        = .TRUE.
USER_COMPXS%ARRAYSDEFINED  = .TRUE.
USER_COMPXS%NGROUP         = NGROUP
USER_COMPXS%NFAM           = NFAM
USER_COMPXS%NCMP           = NCMP
USER_COMPXS%ISCHI          = ISCHI
USER_COMPXS%MAXORD         = MAXORD
USER_COMPXS%LEGTERMS       = LEGTERMS
USER_COMPXS%NKFAM = 0

CALL COMPXS_ZERO(USER_COMPXS)

END SUBROUTINE COMPXS_DEFINE

!---------------------------------------------------------------------------------------------------------------------------------
!  COMPXS_COPY(SOURCE_COMPXS,DESTINATION_COMPXS)
!  Copies the data from SOURCE_COMPXS to DESTINATION_COMPXS
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE COMPXS_COPY(SOURCE_COMPXS,DESTINATION_COMPXS)
IMPLICIT NONE
! Passed in
TYPE (COMPXS_DATA) SOURCE_COMPXS,DESTINATION_COMPXS
! Local
INTEGER I,J,K,L,M

100 FORMAT('[COMPXS]...SORRY, BUT I MUST STOP')
101 FORMAT('[COMPXS]',107('.'))
102 FORMAT('')
105 FORMAT('[COMPXS]...THERE WAS A FATAL ERROR THAT OCCURED IN (ALLOCATE)',54('.'))

IF (SOURCE_COMPXS%DEFINED) THEN
   ! Check to see that it is the correct size
   IF ( DESTINATION_COMPXS%DEFINED .AND. &
        ((DESTINATION_COMPXS%NCMP   .NE. SOURCE_COMPXS%NCMP)  .OR. &
         (DESTINATION_COMPXS%ISCHI  .NE. SOURCE_COMPXS%ISCHI)  .OR. &
         (DESTINATION_COMPXS%MAXORD .NE. SOURCE_COMPXS%MAXORD) .OR. &
         (DESTINATION_COMPXS%NGROUP .NE. SOURCE_COMPXS%NGROUP)      ) ) CALL COMPXS_Void(DESTINATION_COMPXS)
   IF (.NOT. DESTINATION_COMPXS%DEFINED) THEN
      I = SOURCE_COMPXS%NCMP
      J = SOURCE_COMPXS%NGROUP
      K = SOURCE_COMPXS%ISCHI
      L = SOURCE_COMPXS%NFAM
      M = SOURCE_COMPXS%MAXORD
      CALL COMPXS_DEFINE(DESTINATION_COMPXS,I,J,K,L,M)
   END IF
   ! CARD TYPE 1
   DESTINATION_COMPXS%NCMP     = SOURCE_COMPXS%NCMP    
   DESTINATION_COMPXS%NGROUP   = SOURCE_COMPXS%NGROUP  
   DESTINATION_COMPXS%ISCHI    = SOURCE_COMPXS%ISCHI   
   DESTINATION_COMPXS%NFCMP    = SOURCE_COMPXS%NFCMP   
   DESTINATION_COMPXS%MAXUP    = SOURCE_COMPXS%MAXUP   
   DESTINATION_COMPXS%MAXDN    = SOURCE_COMPXS%MAXDN   
   DESTINATION_COMPXS%NFAM     = SOURCE_COMPXS%NFAM    
   DESTINATION_COMPXS%MAXORD   = SOURCE_COMPXS%MAXORD  
   DESTINATION_COMPXS%LEGTERMS = SOURCE_COMPXS%LEGTERMS
   DESTINATION_COMPXS%NDUM2    = SOURCE_COMPXS%NDUM2   
   DESTINATION_COMPXS%NDUM3    = SOURCE_COMPXS%NDUM3   
   ! CARD TYPE 2
   DO I = 1,DESTINATION_COMPXS%NGROUP
      DESTINATION_COMPXS%VEL(I)  = SOURCE_COMPXS%VEL(I)
      DESTINATION_COMPXS%EMAX(I) = SOURCE_COMPXS%EMAX(I)
      DO J = 1,DESTINATION_COMPXS%ISCHI
         DESTINATION_COMPXS%CHI(J,I) = SOURCE_COMPXS%CHI(J,I)
      END DO
   END DO
   I = DESTINATION_COMPXS%NGROUP
   DESTINATION_COMPXS%EMAX(I+1) = SOURCE_COMPXS%EMAX(I+1)
   DO I = 1,DESTINATION_COMPXS%NFAM
      DESTINATION_COMPXS%FLAM(I) = SOURCE_COMPXS%FLAM(I)
      DO J = 1,DESTINATION_COMPXS%NGROUP
         DESTINATION_COMPXS%CHID(J,I) = SOURCE_COMPXS%CHID(J,I)
      END DO
   END DO
   DO L = 1,DESTINATION_COMPXS%NCMP
      DESTINATION_COMPXS%NKFAM(L) = DESTINATION_COMPXS%NKFAM(L)
      ! CARD TYPE 3   COMPOSITION SPECIFICATIONS
      DESTINATION_COMPXS%ICHI(L) = SOURCE_COMPXS%ICHI(L)
      DO I = 1,DESTINATION_COMPXS%NGROUP
         DESTINATION_COMPXS%NUP(I,L)  = SOURCE_COMPXS%NUP(I,L)
         DESTINATION_COMPXS%NDN(I,L)  = SOURCE_COMPXS%NDN(I,L)
         ! CARD TYPE 4   COMPOSITION MACROSCOPIC GROUP CROSS SECTIONS
         DESTINATION_COMPXS%XA(I,L)   = SOURCE_COMPXS%XA(I,L)      
         DESTINATION_COMPXS%XTOT(I,L) = SOURCE_COMPXS%XTOT(I,L)    
         DESTINATION_COMPXS%XREM(I,L) = SOURCE_COMPXS%XREM(I,L)    
         DESTINATION_COMPXS%XTR(I,L)  = SOURCE_COMPXS%XTR(I,L)     
         DESTINATION_COMPXS%XF(I,L)   = SOURCE_COMPXS%XF(I,L)      
         DESTINATION_COMPXS%XNF(I,L)  = SOURCE_COMPXS%XNF(I,L)     
         DESTINATION_COMPXS%PC(I,L)   = SOURCE_COMPXS%PC(I,L)      
         DESTINATION_COMPXS%A1(I,L)   = SOURCE_COMPXS%A1(I,L)      
         DESTINATION_COMPXS%B1(I,L)   = SOURCE_COMPXS%B1(I,L)      
         DESTINATION_COMPXS%A2(I,L)   = SOURCE_COMPXS%A2(I,L)      
         DESTINATION_COMPXS%B2(I,L)   = SOURCE_COMPXS%B2(I,L)      
         DESTINATION_COMPXS%A3(I,L)   = SOURCE_COMPXS%A3(I,L)      
         DESTINATION_COMPXS%B3(I,L)   = SOURCE_COMPXS%B3(I,L)      
         DESTINATION_COMPXS%XN2N(I,L) = SOURCE_COMPXS%XN2N(I,L)
         DO J = 1,DESTINATION_COMPXS%NGROUP !DESTINATION_COMPXS%ICHI(L)
            DESTINATION_COMPXS%ZONE_CHI(I,J,L)  = SOURCE_COMPXS%ZONE_CHI(I,J,L)
         END DO
         DO K = 1,DESTINATION_COMPXS%LEGTERMS
            DO J = 1,DESTINATION_COMPXS%NGROUP
               DESTINATION_COMPXS%XSCAT(J,I,K,L) = SOURCE_COMPXS%XSCAT(J,I,K,L)
            END DO
         END DO
         DO J = 1,DESTINATION_COMPXS%NFAM
            DESTINATION_COMPXS%SNUDEL(J,I,L) = SOURCE_COMPXS%SNUDEL(J,I,L)
         END DO
      END DO
      DO I = 1,DESTINATION_COMPXS%NFAM
         DESTINATION_COMPXS%NUMFAM(I,L) = SOURCE_COMPXS%NUMFAM(I,L)
      END DO
      ! CARD TYPE 5   POWER CONVERSION FACTORS
      DESTINATION_COMPXS%FPWS(L) = SOURCE_COMPXS%FPWS(L)
      DESTINATION_COMPXS%CPWS(L) = SOURCE_COMPXS%CPWS(L)
   END DO
END IF

END SUBROUTINE COMPXS_COPY

!---------------------------------------------------------------------------------------------------------------------------------
!  COMPXS_ADDSCALE(TargetFactor,Target_COMPXS,AddOnFactor,AddOn_COMPXS)
!  Scales the data in Target_COMPXS and adds on AddOn_COMPXS with a factor
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE COMPXS_ADDSCALE(TargetFactor,Target_COMPXS,AddOnFactor,AddOn_COMPXS)
IMPLICIT NONE
! Passed in
REAL*8             AddOnFactor,TargetFactor
TYPE (COMPXS_DATA) AddOn_COMPXS,Target_COMPXS
! Local
INTEGER I,J,K,L,M
LOGICAL Failure

100 FORMAT('[COMPXS]...SORRY, BUT I MUST STOP')
101 FORMAT('[COMPXS]',107('.'))
102 FORMAT('')
105 FORMAT('[COMPXS]...THERE WAS A FATAL ERROR THAT OCCURED IN (ALLOCATE)',54('.'))

Failure = .TRUE.
IF (AddOn_COMPXS%DEFINED .AND. Target_COMPXS%Defined) THEN
   ! Check to see that they are identical sizes
   IF ((Target_COMPXS%NCMP   .NE. AddOn_COMPXS%NCMP)   .OR. &
       (Target_COMPXS%ISCHI  .NE. AddOn_COMPXS%ISCHI)  .OR. &
       (Target_COMPXS%MAXORD .NE. AddOn_COMPXS%MAXORD) .OR. &
       (Target_COMPXS%NGROUP .NE. AddOn_COMPXS%NGROUP)      ) THEN
   ELSE ! Same sizes
      Failure = .FALSE. 
      ! CARD TYPE 2
      DO I = 1,Target_COMPXS%NGROUP
         DO J = 1,Target_COMPXS%ISCHI
            Target_COMPXS%CHI(J,I) = TargetFactor*Target_COMPXS%CHI(J,I) + AddOnFactor*AddOn_COMPXS%CHI(J,I)
         END DO
      END DO
      DO L = 1,Target_COMPXS%NCMP
         ! CARD TYPE 3   COMPOSITION SPECIFICATIONS
         !  ICHI=-1      IF COMPOSITION USES THE SET-WIDE PROMPT CHI GIVEN IN SET CHI RECORD (BELOW). 
         !  ICHI=0       IF COMPOSITION IS NOT FISSIONABLE. 
         !  ICHI=1       FOR COMPOSITION PROMPT CHI VECTOR. 
         !  ICHI=NGROUP FOR COMPOSITION PROMPT CHI MATRIX.
         IF      ((Target_COMPXS%ICHI(L) .EQ. 0) .AND. (AddOn_COMPXS%ICHI(L) .GT.  0)) THEN  ! Target has no fission but add on has isotopic fission
            Target_COMPXS%ICHI(L) = AddOn_COMPXS%ICHI(L)
            DO I = 1,Target_COMPXS%NGROUP
               DO J = 1,Target_COMPXS%NGROUP
                  Target_COMPXS%ZONE_CHI(J,I,L) = AddOn_COMPXS%ZONE_CHI(J,I,L)
               END DO
            END DO
         ELSE IF ((Target_COMPXS%ICHI(L) .EQ. 0) .AND. (AddOn_COMPXS%ICHI(L) .EQ. -1)) THEN  ! Target has no fission but add on used file chi
            Target_COMPXS%ICHI(L) = -1
         ELSE IF ((Target_COMPXS%ICHI(L) .EQ. 1) .AND. (AddOn_COMPXS%ICHI(L) .EQ. -1)) THEN  ! Target has fission but add on used file chi
            DO I = 1,Target_COMPXS%NGROUP
               Target_COMPXS%ZONE_CHI(I,1,L) = TargetFactor*Target_COMPXS%ZONE_CHI(I,1,L) + AddOnFactor*AddOn_COMPXS%CHI(1,I)
            END DO
         ELSE IF ((Target_COMPXS%ICHI(L) .EQ. 1) .AND. (AddOn_COMPXS%ICHI(L) .EQ.  1)) THEN  ! Both have isotopic fission
            DO I = 1,Target_COMPXS%NGROUP
               DO J = 1,Target_COMPXS%NGROUP
                  Target_COMPXS%ZONE_CHI(J,I,L)=TargetFactor*Target_COMPXS%ZONE_CHI(J,I,L)+AddOnFactor*AddOn_COMPXS%ZONE_CHI(J,I,L)
               END DO
            END DO
         END IF
         DO I = 1,Target_COMPXS%NGROUP
            IF (AddOn_COMPXS%NUP(I,L) .GT. Target_COMPXS%NUP(I,L)) Target_COMPXS%NUP(I,L) = AddOn_COMPXS%NUP(I,L)
            IF (AddOn_COMPXS%NDN(I,L) .GT. Target_COMPXS%NDN(I,L)) Target_COMPXS%NDN(I,L) = AddOn_COMPXS%NDN(I,L)
            ! CARD TYPE 4   COMPOSITION MACROSCOPIC GROUP CROSS SECTIONS
            Target_COMPXS%XA(I,L)   =TargetFactor*Target_COMPXS%XA(I,L)   +AddOnFactor*AddOn_COMPXS%XA(I,L)      
            Target_COMPXS%XTOT(I,L) =TargetFactor*Target_COMPXS%XTOT(I,L) +AddOnFactor*AddOn_COMPXS%XTOT(I,L)    
            Target_COMPXS%XREM(I,L) =TargetFactor*Target_COMPXS%XREM(I,L) +AddOnFactor*AddOn_COMPXS%XREM(I,L)    
            Target_COMPXS%XTR(I,L)  =TargetFactor*Target_COMPXS%XTR(I,L)  +AddOnFactor*AddOn_COMPXS%XTR(I,L)     
            Target_COMPXS%XF(I,L)   =TargetFactor*Target_COMPXS%XF(I,L)   +AddOnFactor*AddOn_COMPXS%XF(I,L)      
            Target_COMPXS%XNF(I,L)  =TargetFactor*Target_COMPXS%XNF(I,L)  +AddOnFactor*AddOn_COMPXS%XNF(I,L)     
            Target_COMPXS%PC(I,L)   =TargetFactor*Target_COMPXS%PC(I,L)   +AddOnFactor*AddOn_COMPXS%PC(I,L)      
            Target_COMPXS%A1(I,L)   =TargetFactor*Target_COMPXS%A1(I,L)   +AddOnFactor*AddOn_COMPXS%A1(I,L)      
            Target_COMPXS%B1(I,L)   =TargetFactor*Target_COMPXS%B1(I,L)   +AddOnFactor*AddOn_COMPXS%B1(I,L)      
            Target_COMPXS%A2(I,L)   =TargetFactor*Target_COMPXS%A2(I,L)   +AddOnFactor*AddOn_COMPXS%A2(I,L)      
            Target_COMPXS%B2(I,L)   =TargetFactor*Target_COMPXS%B2(I,L)   +AddOnFactor*AddOn_COMPXS%B2(I,L)      
            Target_COMPXS%A3(I,L)   =TargetFactor*Target_COMPXS%A3(I,L)   +AddOnFactor*AddOn_COMPXS%A3(I,L)      
            Target_COMPXS%B3(I,L)   =TargetFactor*Target_COMPXS%B3(I,L)   +AddOnFactor*AddOn_COMPXS%B3(I,L)      
            Target_COMPXS%XN2N(I,L) =TargetFactor*Target_COMPXS%XN2N(I,L) +AddOnFactor*AddOn_COMPXS%XN2N(I,L)
            DO K = 1,Target_COMPXS%LEGTERMS
               DO J = 1,Target_COMPXS%NGROUP
                  Target_COMPXS%XSCAT(J,I,K,L) = TargetFactor*Target_COMPXS%XSCAT(J,I,K,L) &
                                               + AddOnFactor*AddOn_COMPXS%XSCAT(J,I,K,L)
               END DO
            END DO
         END DO
         ! CARD TYPE 5   POWER CONVERSION FACTORS
         Target_COMPXS%FPWS(L) = TargetFactor*Target_COMPXS%FPWS(L)+AddOnFactor*AddOn_COMPXS%FPWS(L)
         Target_COMPXS%CPWS(L) = TargetFactor*Target_COMPXS%CPWS(L)+AddOnFactor*AddOn_COMPXS%CPWS(L)
      END DO ! L
   END IF 
END IF

END SUBROUTINE COMPXS_ADDSCALE

!---------------------------------------------------------------------------------------------------------------------------------
!  COMPXS_ZERO(USER_COMPXS)
!  Zeros all of the arrays in the COMPXS data structure
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE COMPXS_ZERO(USER_COMPXS)
IMPLICIT NONE
! Passed in
TYPE (COMPXS_DATA) USER_COMPXS                   ! A user data variable to be defined by reading in a COMPXS file
! Local
INTEGER I,J,K,L

100 FORMAT('[COMPXS]...SORRY, BUT I MUST STOP')
101 FORMAT('[COMPXS]',107('.'))
102 FORMAT('')
105 FORMAT('[COMPXS]...THERE WAS A FATAL ERROR THAT OCCURED IN (ALLOCATE)',54('.'))

IF (USER_COMPXS%DEFINED) THEN
   ! CARD TYPE 2
   DO I = 1,USER_COMPXS%NGROUP
      USER_COMPXS%VEL(I)  = 0.0d0 !VEL(I)
      USER_COMPXS%EMAX(I) = 0.0d0 !EMAX(I)
      DO J = 1,USER_COMPXS%ISCHI
         USER_COMPXS%CHI(J,I) = 0.0d0 !CHI(J,I)
      END DO
   END DO
   I = USER_COMPXS%NGROUP
   USER_COMPXS%EMAX(I+1) = 0.0d0 !EMAX(I+1)
   DO I = 1,USER_COMPXS%NFAM
      USER_COMPXS%FLAM(I) = 0.0d0 !FLAM(I)
      DO J = 1,USER_COMPXS%NGROUP
         USER_COMPXS%CHID(J,I) = 0.0d0 !CHID(J,I)
      END DO
   END DO
   DO L = 1,USER_COMPXS%NCMP
      USER_COMPXS%NKFAM(L) = USER_COMPXS%NKFAM(L)
      ! CARD TYPE 3   COMPOSITION SPECIFICATIONS
      USER_COMPXS%ICHI(L) = 0.0d0 !ICHI(L)
      DO I = 1,USER_COMPXS%NGROUP
         USER_COMPXS%NUP(I,L)  = 0.0d0 !NUP(I,L)
         USER_COMPXS%NDN(I,L)  = 0.0d0 !NDN(I,L)
         ! CARD TYPE 4   COMPOSITION MACROSCOPIC GROUP CROSS SECTIONS
         USER_COMPXS%XA(I,L)   = 0.0d0 !XA(I,L)      
         USER_COMPXS%XTOT(I,L) = 0.0d0 !XTOT(I,L)    
         USER_COMPXS%XREM(I,L) = 0.0d0 !XREM(I,L)    
         USER_COMPXS%XTR(I,L)  = 0.0d0 !XTR(I,L)     
         USER_COMPXS%XF(I,L)   = 0.0d0 !XF(I,L)      
         USER_COMPXS%XNF(I,L)  = 0.0d0 !XNF(I,L)     
         USER_COMPXS%PC(I,L)   = 0.0d0 !PC(I,L)      
         USER_COMPXS%A1(I,L)   = 0.0d0 !A1(I,L)      
         USER_COMPXS%B1(I,L)   = 0.0d0 !B1(I,L)      
         USER_COMPXS%A2(I,L)   = 0.0d0 !A2(I,L)      
         USER_COMPXS%B2(I,L)   = 0.0d0 !B2(I,L)      
         USER_COMPXS%A3(I,L)   = 0.0d0 !A3(I,L)      
         USER_COMPXS%B3(I,L)   = 0.0d0 !B3(I,L)      
         USER_COMPXS%XN2N(I,L) = 0.0d0 !XN2N(I,L)
         DO J = 1,USER_COMPXS%NGROUP !USER_COMPXS%ICHI(L)
            USER_COMPXS%ZONE_CHI(I,J,L)  = 0.0d0 !ZONE_CHI(I,J,L)
         END DO
         DO K = 1,USER_COMPXS%LEGTERMS
            DO J = 1,USER_COMPXS%NGROUP
               USER_COMPXS%XSCAT(J,I,K,L) = 0.0d0 !XSCAT(J,I,K,L)
            END DO
         END DO
         DO J = 1,USER_COMPXS%NFAM
            USER_COMPXS%SNUDEL(J,I,L) = 0.0d0 !SNUDEL(J,I,L)
         END DO
      END DO
      DO I = 1,USER_COMPXS%NFAM
         USER_COMPXS%NUMFAM(I,L) = 0.0d0 !NUMFAM(I,L)
      END DO
      ! CARD TYPE 5   POWER CONVERSION FACTORS
      USER_COMPXS%FPWS(L) = 0.0d0 !FPWS(L)
      USER_COMPXS%CPWS(L) = 0.0d0 !CPWS(L)
   END DO
END IF

END SUBROUTINE COMPXS_ZERO

!---------------------------------------------------------------------------------------------------------------------------------
!  COMPXS_VOID(USER_COMPXS)
!  Deallocates the data type
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE COMPXS_VOID(USER_COMPXS)
IMPLICIT NONE
! passed in
TYPE (COMPXS_DATA) USER_COMPXS                   ! A user data variable to be defined by reading in a COMPXS file
! local
DIF3D_Int IOS

100 FORMAT('[COMPXS]...SORRY, BUT I MUST STOP')
101 FORMAT('[COMPXS]',107('.'))
102 FORMAT('')
105 FORMAT('[COMPXS]...THERE WAS A FATAL ERROR THAT OCCURED IN (VOID_DATATYPE)',49('.'))

IF (USER_COMPXS%ARRAYSDEFINED) THEN
   DEALLOCATE( USER_COMPXS%VEL,USER_COMPXS%EMAX,                                    &
               USER_COMPXS%CHID,USER_COMPXS%FLAM,USER_COMPXS%NKFAM,                 &
               USER_COMPXS%NUMFAM,USER_COMPXS%SNUDEL,                               &
               USER_COMPXS%CHI, USER_COMPXS%ICHI,USER_COMPXS%ZONE_CHI,              &
               USER_COMPXS%NUP, USER_COMPXS%NDN,                                    &
               USER_COMPXS%XA,USER_COMPXS%XTOT,USER_COMPXS%XREM,USER_COMPXS%XTR,    &
               USER_COMPXS%XF,USER_COMPXS%XNF,USER_COMPXS%XSCAT,                    &
               USER_COMPXS%XN2N,USER_COMPXS%PC,                                     &
               USER_COMPXS%A1,USER_COMPXS%A2,USER_COMPXS%A3,                        &
               USER_COMPXS%B1,USER_COMPXS%B2,USER_COMPXS%B3,                        &
               USER_COMPXS%FPWS,USER_COMPXS%CPWS,                                   &
               STAT = IOS                                                           )
   IF (IOS .NE. 0) THEN
      WRITE(MODULE_OUT,105)
      WRITE(MODULE_OUT,'("[COMPXS]...FAILED TO DEALLOCATE COMPXS DATA",72("."))')
      WRITE(MODULE_OUT,100)
      CALL Abort
   END IF
   USER_COMPXS%ARRAYSDEFINED  = .FALSE.
END IF

USER_COMPXS%DEFINED  = .FALSE.
USER_COMPXS%NCMP     = 0
USER_COMPXS%NGROUP   = 0
USER_COMPXS%ISCHI    = 0
USER_COMPXS%NFCMP    = 0
USER_COMPXS%MAXUP    = 0
USER_COMPXS%MAXDN    = 0
USER_COMPXS%NFAM     = 0
USER_COMPXS%MAXORD   = 0
USER_COMPXS%LEGTERMS = 0
USER_COMPXS%NDUM2    = 0
USER_COMPXS%NDUM3    = 0

END SUBROUTINE COMPXS_VOID

!---------------------------------------------------------------------------------------------------------------------------------
!  COMPXS_UPDATEPROPERTIES(USER_COMPXS)
!  Update several property variable pieces of the existing COMPXS data set before exporting it
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE COMPXS_UPDATEPROPERTIES(USER_COMPXS)
IMPLICIT NONE
! Passed in
TYPE (COMPXS_DATA) USER_COMPXS                   ! A user data variable to be defined by reading in a COMPXS file
! Local
DIF3D_Int CURC,G,GG,J,I
DIF3D_Int NFCMP,MAXORD
DIF3D_Log FOUNDFISSION,FOUNDCHI,FOUNDCHIMATRIX

! This routine will determine the following pieces of information
!   DIF3D_Int NFCMP                                               ! NUMBER OF FISSIONABLE COMPOSITIONS.
!   DIF3D_Int MAXUP                                               ! MAXIMUM NUMBER OF GROUPS OF UPSCATTERING FOR THE SET.
!   DIF3D_Int MAXDN                                               ! MAXIMUM NUMBER OF GROUPS OF DOWNSCATTERING FOR THE SET.
!   DIF3D_Int MAXORD                                              ! ANISOTROPIC SCATTERING ORDER
!   DIF3D_Int LEGTERMS                                            ! NUMBER OF LEGENDRE TERMS ! MAXORD + 1
!   DIF3D_Int,          DIMENSION(:),      POINTER :: ICHI    ! NCMP         ; PROMPT FISSION SPECTRUM FLAG FOR THIS COMPOSITION. 
!   !                                                                          ; ICHI=-1      IF COMPOSITION USES THE SET-WIDE PROMPT CHI GIVEN IN SET CHI RECORD (BELOW). 
!   !                                                                          ; ICHI=0       IF COMPOSITION IS NOT FISSIONABLE. 
!   !                                                                          ; ICHI=1       FOR COMPOSITION PROMPT CHI VECTOR. 
!   !                                                                          ; ICHI=NGROUP FOR COMPOSITION PROMPT CHI MATRIX.
!   DIF3D_Int,          DIMENSION(:,:),    POINTER :: NUP     ! NGROUP,NCMP  ; NUMBER OF GROUPS OF UPSCATTERING INTO GROUP I FROM LOWER ENERGY GROUPS FOR THE CURRENT COMPOSITION
!   DIF3D_Int,          DIMENSION(:,:),    POINTER :: NDN     ! NGROUP,NCMP  ; NUMBER OF GROUPS OF DOWNSCATTERING INTO GROUP I FROM HIGHER ENERGY GROUPS FOR THE CURRENT COMPOSITION
! Zero out the pieces that will be determined from the existing data
NFCMP = 0
USER_COMPXS%MAXUP = 0
USER_COMPXS%MAXDN = 0
MAXORD = 0
USER_COMPXS%ICHI = 0
USER_COMPXS%NDN = 0
USER_COMPXS%NUP = 0

DO CURC = 1,USER_COMPXS%NCMP
   ! Sweep through compute the number of fissionable compositions and set ICHI if non-zero XF and the subsequent presence of ZONE_CHI
   FOUNDFISSION = .FALSE.
   JUNK1: DO G = 1,USER_COMPXS%NGROUP
      IF (USER_COMPXS%XF(G,CURC) .NE. 0.0D0) THEN
         FOUNDFISSION = .TRUE.
         EXIT JUNK1
      END IF
   END DO JUNK1
   IF (FOUNDFISSION) THEN
      NFCMP = NFCMP + 1
      FOUNDCHI = .FALSE.
      FOUNDCHIMATRIX = .FALSE.
      JUNK2: DO G = 1,USER_COMPXS%NGROUP
         IF (USER_COMPXS%ZONE_CHI(G,1,CURC) .NE. 0.0D0) FOUNDCHI = .TRUE.  ! At least chi exists
         DO GG = 2,USER_COMPXS%NGROUP
            IF (USER_COMPXS%ZONE_CHI(G,GG,CURC) .NE. 0.0D0) THEN
               FOUNDCHIMATRIX = .TRUE.    ! Entire chi matrix exists
               EXIT JUNK2
            END IF
         END DO
      END DO JUNK2

      IF (FOUNDCHIMATRIX) THEN
         USER_COMPXS%ICHI(CURC) = USER_COMPXS%NGROUP  ! Local Chi matrix exists
      ELSE IF (FOUNDCHI) THEN
         USER_COMPXS%ICHI(CURC) = 1    ! Only the Chi vector is present
      ELSE
         USER_COMPXS%ICHI(CURC) = -1   ! Use the set Chi
      END IF
   ELSE
      USER_COMPXS%ICHI(CURC) = 0 ! Not fissionable
   END IF
   ! MAXORD is only checked if it has a chance for being changed and within group scattering must be non-zero for it to be changed
   IF (MAXORD .LT. USER_COMPXS%MAXORD) THEN
      JUNK3: DO I = USER_COMPXS%MAXORD,MAXORD+1,-1    ! 2,1 as an example
         DO G = 1,USER_COMPXS%NGROUP
            IF (USER_COMPXS%XSCAT(G,G,I+1,CURC) .NE. 0.0D0) THEN
               MAXORD = I  ! This is clearly the new highest order
               EXIT JUNK3
            END IF
         END DO
      END DO JUNK3
   END IF
   ! Sweep through all compositions and compute new scattering bandwidths and maximum bandwidths based upon isotropic scattering
   DO G = 1,USER_COMPXS%NGROUP   ! SCATTERING TO GROUP G  SCAT(GG->G)
      USER_COMPXS%NDN(G,CURC) = 0
      GG = 0
      DO ! GG = 1,G-1                                         ! DOWN FROM GROUP GG
         GG = GG + 1
         ! WRITE(MODULE_OUT,'("SCAT(",I5,",",I5,") =",E13.5)') GG,G,USER_COMPXS%XSCAT(GG,G,1,CURC)
         IF (GG .EQ. G) EXIT                                  ! hit limiting condition on GG so exit
         IF (USER_COMPXS%XSCAT(GG,G,1,CURC) .NE. 0.0D0) EXIT  ! Found first non-zero number so exit
      END DO
      ! WRITE(MODULE_OUT,'("GG = ",I5," G-GG=",I5)')GG,G-GG
      USER_COMPXS%NDN(G,CURC) = G - GG    ! Length of scattering block (cannot be less than zero
      IF (USER_COMPXS%NDN(G,CURC) .GT. USER_COMPXS%MAXDN) USER_COMPXS%MAXDN = USER_COMPXS%NDN(G,CURC)

      USER_COMPXS%NUP(G,CURC) = 0
      GG = 0
      DO ! GG = 1,USER_COMPXS%NGROUP-G                        ! UP FROM GROUP GG
         GG = GG + 1
         J = USER_COMPXS%NGROUP+1-GG                          ! group location
         ! WRITE(MODULE_OUT,'("SCAT(",I5,",",I5,") =",E13.5)') J,G,USER_COMPXS%XSCAT(J,G,1,CURC)
         IF (J .EQ. G) EXIT                                   ! hit limiting condition on GG so exit
         IF (USER_COMPXS%XSCAT(J,G,1,CURC) .NE. 0.0D0) EXIT   ! Found first non-zero number so exit
      END DO
      ! WRITE(MODULE_OUT,'("GG = ",I5," J = ",I5)')GG,J
      USER_COMPXS%NUP(G,CURC) = J-G    ! Length of scattering block (cannot be less than zero
      IF (USER_COMPXS%NUP(G,CURC) .GT. USER_COMPXS%MAXUP) USER_COMPXS%MAXUP = USER_COMPXS%NUP(G,CURC)
   END DO
END DO

! Set the variable data
USER_COMPXS%NFCMP    = NFCMP
IF (MAXORD .GT. USER_COMPXS%MAXORD) THEN
   USER_COMPXS%MAXORD = MAXORD
   USER_COMPXS%LEGTERMS = MAXORD + 1
END IF

END SUBROUTINE COMPXS_UPDATEPROPERTIES

END MODULE COMPXS_io

