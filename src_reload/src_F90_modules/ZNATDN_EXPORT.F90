!---------------------------------------------------------------------------------------------------------------------------------
!  ZNATDN_EXPORT(USER_ZNATDN,USER_FILENAME)
!  Exports the ZNATDN data file
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ZNATDN_EXPORT(Output_Unit,USER_ZNATDN,USER_FILENAME)
USE ZNATDN_io
USE NE_Kind
#include "DIF3D_Types.h"
IMPLICIT NONE
! Passed in
DIF3D_Int Output_Unit
TYPE (ZNATDN_DATA) USER_ZNATDN  ! A user data variable to be defined by reading in a ZNATDN file
CHARACTER(*) USER_FILENAME
! Local
DIF3D_Int OUTPUTUNIT
DIF3D_Int IOS,I,J,M,N

100 FORMAT('[ZNATDN]...SORRY, BUT I MUST STOP')
101 FORMAT('[ZNATDN]',107('.'))
102 FORMAT('')
105 FORMAT('[ZNATDN]...THERE WAS A FATAL ERROR THAT OCCURED IN (EXPORT)',56('.'))

IF (.NOT. USER_ZNATDN%DEFINED) THEN  ! Data to write doesn't exist
   WRITE(Output_Unit,105)
   WRITE(Output_Unit,'("[ZNATDN]...USER DATA WAS UNDEFINED FOR EXPORT <",A60,">",7("."))') USER_FILENAME
   WRITE(Output_Unit,100)
   CALL Abort
END IF

USER_ZNATDN%FILENAME = ADJUSTL(USER_FILENAME)
OUTPUTUNIT = NE_Kind_GetFreeLogicalUnit()

OPEN(UNIT=OUTPUTUNIT,IOSTAT=IOS,FILE=USER_FILENAME,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
IF (IOS .NE. 0) THEN
   WRITE(Output_Unit,105)
   WRITE(Output_Unit,'("[ZNATDN]...CANNOT OPEN FILE <",A80,">",5("."))') USER_FILENAME
   WRITE(Output_Unit,100)
   CALL Abort
END IF

#ifdef DIF3D_Debug
   WRITE(Output_Unit,101)
   WRITE(Output_Unit,'("[DLAYXS]...DEBUG PRINT OF ZNATDN EXPORT ROUTINE",52("."),A16)') USER_FILENAME
   WRITE(Output_Unit,101)
#endif

! Eliminate blocking
USER_ZNATDN%NBLKAD = 1

! WRITE Card Type 0
WRITE(OUTPUTUNIT) USER_ZNATDN%HNAME,(USER_ZNATDN%HUSE(I),I=1,2),USER_ZNATDN%IVERS
! WRITE Card Type 1
WRITE(OUTPUTUNIT) USER_ZNATDN%TIME,USER_ZNATDN%NCY,USER_ZNATDN%NTZSZ,USER_ZNATDN%NNS,USER_ZNATDN%NBLKAD

#ifdef DIF3D_Debug
   WRITE(Output_Unit,'("[ZNATDN]...HNAME................................",51("."),A16)') USER_ZNATDN%HNAME
   WRITE(Output_Unit,'("[ZNATDN]...HUSE(1)..............................",51("."),A16)') USER_ZNATDN%HUSE(1)
   WRITE(Output_Unit,'("[ZNATDN]...HUSE(2)..............................",51("."),A16)') USER_ZNATDN%HUSE(2)
   WRITE(Output_Unit,'("[ZNATDN]...IVERS................................",51("."),I16)') USER_ZNATDN%IVERS
   WRITE(Output_Unit,'("[ZNATDN]...REFERENCE REAL TIME, DAYS............",51("."),1PE16.9)') USER_ZNATDN%TIME
   WRITE(Output_Unit,'("[ZNATDN]...REFERENCE CYCLE NUMBER...............",51("."),I16)') USER_ZNATDN%NCY
   WRITE(Output_Unit,'("[ZNATDN]...NUMBER OF ZONES PLUS SUBZONES........",51("."),I16)') USER_ZNATDN%NTZSZ
   WRITE(Output_Unit,'("[ZNATDN]...MAXIMUM NUMBER OF NUCLIDES IN ANY SET",51("."),I16)') USER_ZNATDN%NNS
   WRITE(Output_Unit,'("[ZNATDN]...NUMBER OF BLOCKS OF ATOM DENSITY DATA",51("."),I16)') USER_ZNATDN%NBLKAD
#endif

! WRITE Card Type 2
WRITE(OUTPUTUNIT) ((USER_ZNATDN%ADEN(N,J),N=1,USER_ZNATDN%NNS),J=1,USER_ZNATDN%NTZSZ)

#ifdef DIF3D_Debug
   DO J = 1,USER_ZNATDN%NTZSZ
      WRITE(Output_Unit,'("[ZNATDN]...ZONE ATOM DENSITIES FOR MIXTURE",57("."),I16)') J
      WRITE(Output_Unit,'(("[ZNATDN]...",4(I3," [",1PE13.6,"]",3X),16(".")))') (N,USER_ZNATDN%ADEN(N,J),N=1,USER_ZNATDN%NNS)
   END DO
#endif

CLOSE(UNIT=OUTPUTUNIT)

CALL NE_Kind_FREELOGICALUNIT(OUTPUTUNIT)

END SUBROUTINE ZNATDN_EXPORT

