!  This subroutine simplifies the NDXSRF (and thus ZNATDN and LABELS) data sets to consist of a single isotope set and eliminate sub-zones
!  Copyright(c) 2018 Argonne National Laboratory
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ARC_Simplify_NDXSRF(Output_Unit,InNDXSRF,InZNATDN,InLABELS)
#include "DIF3D_Types.h"
!#define Local_Debug
USE ZNATDN_io
USE NDXSRF_io
USE LABELS_io
IMPLICIT NONE
! Passed in
DIF3D_Int                     Output_Unit
TYPE (NDXSRF_Data)            InNDXSRF          ! This structure will be updated with the new isotopes
TYPE (ZNATDN_Data)            InZNATDN          ! This structure will be updated with atom densities equivalent to those of 
TYPE (LABELS_Data)            InLABELS
! Local
TYPE (NDXSRF_Data)            CopyNDXSRF
TYPE (ZNATDN_Data)            CopyZNATDN
TYPE (LABELS_Data)            CopyLABELS
! Local junk
DIF3D_Int I,J,K
DIF3D_Int izone,isubzone,isotope_ISOTXS
DIF3D_Int NON,NSN,NNS,NAN,NZONE,NSZ,IOS
DIF3D_Real SomeReal

100 FORMAT('[ARC]...Sorry, but I must stop')
105 FORMAT('[ARC]...There was a fatal error that occured in (ARC_Simplify_NDXSRF)')
106 FORMAT('[ARC]...There was a non-fatal error that occured in (ARC_Simplify_NDXSRF)')
200 FORMAT('[ARC]...',A80)

IF (.NOT. InNDXSRF%Defined) THEN
   WRITE(Output_Unit,'("[ARC]...NDXSRF data structure not defined")')
   CALL Basic_Abort
END IF

IF (.NOT. InZNATDN%Defined) THEN
   WRITE(Output_Unit,'("[ARC]...ZNATDN data structure not defined")')
   CALL Basic_Abort
END IF

IF (.NOT. InLABELS%Defined) THEN
   WRITE(Output_Unit,'("[ARC]...ZNATDN data structure not defined")')
   CALL Basic_Abort
END IF

#ifdef Local_Debug
   WRITE(Output_Unit,'(//,"[ARC]...Data structures before modifications")')
   CALL NDXSRF_Print(InNDXSRF)
   CALL ZNATDN_Print(InZNATDN)
   CALL LABELS_Print(InLABELS)
   WRITE(Output_Unit,'(//)')
#endif

! The change to labels is very easy
IF (InNDXSRF%NSZ .NE. 0) THEN
   ALLOCATE(CopyLABELS%CMPNAM(InLABELS%NTZSZ),STAT=IOS)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,'("[ARC]...Failed to allocate scratch LABELS memory in ARC_Simplify_NDXSRF")')
      Call Abort
   END IF
   DO I = 1,InNDXSRF%NZONE
      CopyLABELS%CMPNAM(I) = InLABELS%CMPNAM(I)
   END DO
   DEALLOCATE(InLABELS%CMPNAM,STAT=IOS)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,'("[ARC]...Failed to deallocate LABELS memory in ARC_Simplify_NDXSRF")')
      Call Abort
   END IF
   ALLOCATE(InLABELS%CMPNAM(InNDXSRF%NZONE),STAT=IOS)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,'("[ARC]...Failed to allocate LABELS memory in ARC_Simplify_NDXSRF")')
      Call Abort
   END IF
   DO I = 1,InNDXSRF%NZONE
      InLABELS%CMPNAM(I) = CopyLABELS%CMPNAM(I)
   END DO
   DEALLOCATE(CopyLABELS%CMPNAM,STAT=IOS)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,'("[ARC]...Failed to deallocate scratch LABELS memory in ARC_Simplify_NDXSRF")')
      Call Abort
   END IF
END IF ! NSZ <> 0

! Make a copy of the NDXSRF and ZNATDN data
CALL NDXSRF_Copy(InNDXSRF,CopyNDXSRF)
CALL ZNATDN_Copy(InZNATDN,CopyZNATDN)

! Rebuild the NDXSRF and ZNATDN data structures
NON   = InNDXSRF%NON
NSN   = 1 ! InNDXSRF%NSN ! we only want a single isotope set that includes all isotopes in ISOTXS
NNS   = InNDXSRF%NON ! InNDXSRF%NNS  maximum nuclides per set
NAN   = InNDXSRF%NON !               number of different nuclides in the set (some idiot says we needed two definitions of NON and another mapping array)
NZONE = InNDXSRF%NZONE
NSZ   = 0 ! InNDXSRF%NSZ
CALL NDXSRF_DEFINE(InNDXSRF,NON,NSN,NNS,NAN,NZONE,NSZ)
CALL ZNATDN_DEFINE(InZNATDN,NNS,NZONE)

! Copy the base NDXSRF data
DO I = 1,CopyNDXSRF%NON ! old NISO
   InNDXSRF%HNNAME(I) = CopyNDXSRF%HNNAME(I)
   InNDXSRF%HANAME(I) = CopyNDXSRF%HANAME(I)
   InNDXSRF%WPF(I)    = CopyNDXSRF%WPF(I)
   InNDXSRF%ATWT(I)   = CopyNDXSRF%ATWT(I)
   InNDXSRF%NCLN(I)   = CopyNDXSRF%NCLN(I)
END DO

! Define the new NOS and NOR arrays
InNDXSRF%NDXS(1,1) = NON ! Number of isotopes in set
InNDXSRF%NDXS(2,1) = CopyNDXSRF%NDXS(2,1) ! garbage for the rest
InNDXSRF%NDXS(3,1) = CopyNDXSRF%NDXS(3,1)
InNDXSRF%NDXS(4,1) = CopyNDXSRF%NDXS(4,1)
DO J = 1,NON
   InNDXSRF%NOR(J,1) = J ! The ISOTXS index id for isotope J in index set 1
   InNDXSRF%NOS(J,1) = J ! Given ISOTXS id J, what is the order number in this set
END DO

! Copy the ZNATDN data over eliminating all sub-zones
DO izone = 1,CopyNDXSRF%NZONE
   InNDXSRF%VOLZ(izone) = CopyNDXSRF%VOLZ(izone)
   InNDXSRF%VFPA(izone) = 1.0d0 ! Volume fraction
   InNDXSRF%NSPA(izone) = 1     ! isotope set assignment
   K = CopyNDXSRF%NSPA(izone)   ! The index set for the old data
   IF (K .NE. 0) THEN 
      DO I = 1,CopyZNATDN%NNS ! all isotopes in set
         isotope_ISOTXS = CopyNDXSRF%NOS(I,K) ! The ISOTXS isotope associated with the old data
         InZNATDN%ADEN(isotope_ISOTXS,izone) = InZNATDN%ADEN(isotope_ISOTXS,izone) + CopyZNATDN%ADEN(I,izone)*CopyNDXSRF%VFPA(izone)
      END DO
   END IF
   DO isubzone = 1,CopyNDXSRF%NSZ ! subzones
      IF (CopyNDXSRF%NZSZ(isubzone) .EQ. izone) THEN ! sub-zone is in zone
         SomeReal = CopyNDXSRF%VLSA(isubzone)/CopyNDXSRF%VOLZ(izone)
         K = CopyNDXSRF%NSSA(isubzone)   ! The index set for the old data
         IF (K .NE. 0) THEN
            DO I = 1,CopyZNATDN%NNS ! all isotopes in set
               isotope_ISOTXS = CopyNDXSRF%NOS(I,K) ! The ISOTXS isotope associated with the old data
               InZNATDN%ADEN(isotope_ISOTXS,izone) = InZNATDN%ADEN(isotope_ISOTXS,izone)+CopyZNATDN%ADEN(I,NZONE+isubzone)*SomeReal
            END DO
         END IF
      END IF
   END DO
END DO

CALL NDXSRF_Void(CopyNDXSRF)
CALL ZNATDN_Void(CopyZNATDN)

#ifdef Local_Debug
   WRITE(Output_Unit,'(//,"[ARC]...Data structures after modifications")')
   CALL NDXSRF_Print(InNDXSRF)
   CALL ZNATDN_Print(InZNATDN)
   WRITE(Output_Unit,'(//)')
   STOP
#endif

#ifdef Local_Debug
! Print out the compositions
DO izone = 1,InNDXSRF%NZONE
   WRITE(Output_Unit,*)'ZONE ',izone
   DO isotope_ISOTXS = 1,InNDXSRF%NON ! all isotopes in set
      IF (InZNATDN%ADEN(isotope_ISOTXS,izone) .NE. 0.0d0) THEN
         WRITE(Output_Unit,300) InNDXSRF%HNNAME(isotope_ISOTXS),InZNATDN%ADEN(isotope_ISOTXS,izone)
         300 FORMAT('Isotope name ',A8,' density ',1PE13.6)
      END IF
   END DO
END DO
#endif

END SUBROUTINE ARC_Simplify_NDXSRF
