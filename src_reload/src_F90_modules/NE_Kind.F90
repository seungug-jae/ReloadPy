!  Copyright(c) 2005 Argonne National Laboratory
!---------------------------------------------------------------------------------------------------------------------------------
!  A mishmash of stuff used in DIF3D that has no other place to go
!---------------------------------------------------------------------------------------------------------------------------------
!  Subroutines and functions:
!---------------------------------------------------------------------------------------------------------------------------------
!  NE_Kind_SetUnit       !  Sets the output unit
!  NE_Kind_GetFreeLogicalUnit    !  Returns a DIF3D_Log unit number that is currently not connected to a open file
!  NE_Kind_FREELOGICALUNIT       !  Frees a DIF3D_Log unit number that was previously connected to a open file
!  NE_Kind_StandardCharacterFix  !  Takes a character string, left justifies it and uppercases all of the letters, and trims it
!  NE_Kind_DIF3DCharacterFix     !  Takes a character string, left justifies it and removes all spaces and trims it
!  NE_Kind_UpperCaseString       !  Takes a character string and uppercases all of the letters
!  NE_Kind_AddUnderscores        !  Takes a character string and replaces the spaces with underscores '_'
!  NE_Kind_YesOrNo               !  Simple function that returns a character "yes" or "no"
!  NE_Kind_IntToLogical          !  Simple function that converts 0 to false, and 1 to true and returns it in a standard DIF3D_Log
!  NE_Kind_LogicalToInt          !  Simple function that converts false to 0 and true to 1 and returns a standard DIF3D_Int
!  NE_Kind_IntegerAlphabet       !  Takes an integer and converts it to a 6 letter string using the NE_Kind alphabet
!  NE_Kind_AlphabetInteger       !  Takes an 6 letter NE_Kind alphabet word and converts it to an integer using 
!--------------------------------
!  NE_Kind_SetUnit(OUTPUT_UNIT)
!  DIF3D_Int FUNCTION NE_Kind_GetFreeLogicalUnit()
!  NE_Kind_FREELOGICALUNIT(I)
!  NE_Kind_StandardCharacterFix(User_String)
!  NE_Kind_DIF3DCharacterFix(User_String)
!  NE_Kind_UpperCaseString(User_String)
!  NE_Kind_AddUnderscores(User_String)
!  DIF3D_Char_Size FUNCTION NE_Kind_YesOrNo(User_DIF3D_Log)
!  DIF3D_Log FUNCTION NE_Kind_IntToLogical(User_Integer)
!  DIF3D_Int FUNCTION NE_Kind_LogicalToInt(User_DIF3D_Log)
!  NE_Kind_IntegerAlphabet(MyInteger,MyString)
!  NE_Kind_AlphabetInteger(MyInteger,MyString)
!---------------------------------------------------------------------------------------------------------------------------------

MODULE NE_Kind
IMPLICIT NONE
#include "DIF3D_Types.h"
!#include "DIF3D_Macros.h"

DIF3D_Int :: MODULE_OUT = 6   ! Default output unit to use
DIF3D_Int :: NE_Kind_SeekFileUnit = 500    ! Starting unit number to use for seeking free units. I used 300 to put it above DIF3D/REBUS stuff

DIF3D_Int    :: NE_Kind_SizeAlphabet = 63
CHARACTER*63 :: NE_Kind_MyAlphabet = '_0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'


PRIVATE MODULE_OUT,NE_Kind_SeekFileUnit

CONTAINS

!---------------------------------------------------------------------------------------------------------------------------------
!  NE_Kind_SetUnit(OUTPUT_UNIT)
!  Sets the output unit
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE NE_Kind_SetUnit(OUTPUT_UNIT)
DIF3D_Int OUTPUT_UNIT
MODULE_OUT = OUTPUT_UNIT
END SUBROUTINE NE_Kind_SetUnit

!---------------------------------------------------------------------------------------------------------------------------------
!  DIF3D_Int FUNCTION NE_Kind_GetFreeLogicalUnit()
!  Returns a DIF3D_Log unit number that is currently not connected to a open file
!---------------------------------------------------------------------------------------------------------------------------------
DIF3D_Int FUNCTION NE_Kind_GetFreeLogicalUnit()
IMPLICIT NONE
! Local variables
DIF3D_Int IOS
DIF3D_Log ISOPEN
DIF3D_Log L_Attempt
! Function variable
!DIF3D_Int NE_Kind_GetFreeLogicalUnit

100 FORMAT('[NE_Kind]...SORRY, BUT I MUST STOP')

L_Attempt = .FALSE.
! Find an available DIF3D_Log unit number
1000 CONTINUE
DO 
   NE_Kind_SeekFileUnit = NE_Kind_SeekFileUnit + 1
   INQUIRE(UNIT=NE_Kind_SeekFileUnit,OPENED=ISOPEN,IOSTAT=IOS)
   !WRITE(6,*)'I=',NE_Kind_SeekFileUnit,' ISOPEN=',ISOPEN
   IF ((IOS .EQ. 0) .AND. (.NOT. ISOPEN)) EXIT
   IF (NE_Kind_SeekFileUnit .GT. DIF3D_MaxFileUnits) EXIT
END DO

IF ((NE_Kind_SeekFileUnit .GT. DIF3D_MaxFileUnits) .AND. (.NOT. L_Attempt)) THEN
   NE_Kind_SeekFileUnit = 500 ! Restart at the bottom
   L_Attempt = .TRUE. ! Indicate that we cannot look again
   GOTO 1000
END IF

IF (NE_Kind_SeekFileUnit .GT. DIF3D_MaxFileUnits) THEN
   WRITE(MODULE_OUT,'("[NE_Kind]...TOO MANY UNIT NUMBERS (",I6,">",I6,") TRY FREEING SOME USING NE_Kind_FREELOGICALUNIT")') &
         NE_Kind_SeekFileUnit,DIF3D_MaxFileUnits
   WRITE(MODULE_OUT,100)
   CALL Abort
!ELSE
!   WRITE(6,*)' Get SeekFileUnit=',NE_Kind_SeekFileUnit
END IF

NE_Kind_GetFreeLogicalUnit = NE_Kind_SeekFileUnit

END FUNCTION NE_Kind_GetFreeLogicalUnit

!---------------------------------------------------------------------------------------------------------------------------------
!  NE_Kind_FREELOGICALUNIT(I)
!  Frees a DIF3D_Log unit number that was previously connected to a open file
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE NE_Kind_FREELOGICALUNIT(I)
IMPLICIT NONE
! Passed in variables
DIF3D_Int I
! Local variables
DIF3D_Int IOS
DIF3D_Log ISOPEN

100 FORMAT('[NE_Kind]...SORRY, BUT I MUST STOP')

! Compare the freed unit number with the stored number
IF ((I .GT. 6) .AND. (I .LE. DIF3D_MaxFileUnits)) THEN
   IF (I .LE. NE_Kind_SeekFileUnit) NE_Kind_SeekFileUnit = I-1
END IF

IF (NE_Kind_SeekFileUnit .LT. 6) NE_Kind_SeekFileUnit = 6
!WRITE(6,*)'I=',I,' Free SeekFileUnit=',NE_Kind_SeekFileUnit

END SUBROUTINE NE_Kind_FREELOGICALUNIT

!---------------------------------------------------------------------------------------------------------------------------------
!  NE_Kind_StandardCharacterFix(User_String)
!  Takes a character string, left justifies it and uppercases all of the letters, and trims it
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE NE_Kind_StandardCharacterFix(User_String)
IMPLICIT NONE
! Passed in variables
CHARACTER(*) User_String

! If you get a segmentation fault from using this routine then you are passing in a parameter statement dumbass
User_String = ADJUSTL(User_String)                ! Left justify
CALL NE_Kind_UpperCaseString(User_String)    ! uppercase
User_String = TRIM(User_String)                   ! Trim

END SUBROUTINE NE_Kind_StandardCharacterFix

!---------------------------------------------------------------------------------------------------------------------------------
!  NE_Kind_DIF3DCharacterFix(User_String)
!  Takes a character string, left justifies it and removes all spaces and trims it
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE NE_Kind_DIF3DCharacterFix(User_String)
IMPLICIT NONE
! Passed in variables
CHARACTER(*) User_String
DIF3D_Int LENGTH,I,J

! If you get a segmentation fault from using this routine then you are passing in a parameter statement dumbass
User_String = ADJUSTL(User_String)           ! Left justify

LENGTH = LEN(User_String)
IF (LENGTH .GT. 0) THEN
   J = 0
   DO I = 1,LENGTH
      IF (User_String(I:I) .NE. ' ') THEN
         J = J + 1
         User_String(J:J) = User_String(I:I)
      END IF
   END DO
   DO I = J+1,LENGTH
      User_String(I:I) = ''
   END DO
END IF
User_String = TRIM(User_String)              ! Trim

END SUBROUTINE NE_Kind_DIF3DCharacterFix

!---------------------------------------------------------------------------------------------------------------------------------
!  NE_Kind_UpperCaseString(User_String)
!  Takes a character string and uppercases all of the letters
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE NE_Kind_UpperCaseString(User_String)
IMPLICIT NONE
! Passed in variables
CHARACTER(*) User_String
! Local variables
DIF3D_Int LENGTH,I,J

LENGTH = LEN(User_String)

IF (LENGTH .GT. 0) THEN
   DO I = 1,LENGTH
      J = IACHAR(User_String(I:I))
      IF ((J .GE. 97) .AND. (J .LE. 122)) User_String(I:I) = ACHAR(J-32)
   END DO
END IF

END SUBROUTINE NE_Kind_UpperCaseString

!---------------------------------------------------------------------------------------------------------------------------------
!  NE_Kind_AddUnderscores(User_String)
!  Takes a character string and replaces the spaces with underscores '_'
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE NE_Kind_AddUnderscores(User_String)
IMPLICIT NONE
! Passed in variables
CHARACTER(*) User_String
! Local variables
DIF3D_Int LENGTH,I,J

LENGTH = LEN(User_String)

IF (LENGTH .GT. 0) THEN
   ! Find the last non-blank character
   Search: DO J = LENGTH,1,-1
      IF (User_String(J:J) .NE. ' ') EXIT Search
   END DO Search

   DO I = 1,J
      IF (User_String(I:I) .EQ. ' ') User_String(I:I) = '_'
   END DO
END IF

END SUBROUTINE NE_Kind_AddUnderscores

!---------------------------------------------------------------------------------------------------------------------------------
!  DIF3D_Char_Size FUNCTION NE_Kind_YesOrNo(User_DIF3D_Log)
!  Simple function that returns a character "yes" or "no"
!---------------------------------------------------------------------------------------------------------------------------------
DIF3D_Char FUNCTION NE_Kind_YesOrNo(User_DIF3D_Log)
IMPLICIT NONE
DIF3D_Log :: User_DIF3D_Log    ! Passed in
IF (User_DIF3D_Log) THEN
   NE_Kind_YesOrNo = DIF3D_Definition_Yes;
ELSE
   NE_Kind_YesOrNo = DIF3D_Definition_No;
END IF
END FUNCTION NE_Kind_YesOrNo

!---------------------------------------------------------------------------------------------------------------------------------
!  DIF3D_Log FUNCTION NE_Kind_IntToLogical(User_Integer)
!  Simple function that converts 0 to false, and 1 to true and returns it in a standard DIF3D_Log
!---------------------------------------------------------------------------------------------------------------------------------
DIF3D_Log FUNCTION NE_Kind_IntToLogical(User_Integer)
IMPLICIT NONE
DIF3D_Int :: User_Integer    ! Passed in
IF (User_Integer .EQ. 0) THEN
   NE_Kind_IntToLogical = .FALSE.;
ELSE
   NE_Kind_IntToLogical = .TRUE.;
END IF
END FUNCTION NE_Kind_IntToLogical

!---------------------------------------------------------------------------------------------------------------------------------
!  DIF3D_Int FUNCTION NE_Kind_LogicalToInt(User_DIF3D_Log)
!  Simple function that converts false to 0 and true to 1 and returns a standard DIF3D_Int
!---------------------------------------------------------------------------------------------------------------------------------
DIF3D_Int FUNCTION NE_Kind_LogicalToInt(User_DIF3D_Log)
IMPLICIT NONE
DIF3D_Log :: User_DIF3D_Log    ! Passed in
IF (User_DIF3D_Log) THEN
   NE_Kind_LogicalToInt = 1;
ELSE
   NE_Kind_LogicalToInt = 0;
END IF
END FUNCTION NE_Kind_LogicalToInt

!---------------------------------------------------------------------------------------------------------------------------------
!  NE_Kind_IntegerAlphabet(MyInteger,MyString)
!  Takes an integer and converts it to a 6 letter string using the NE_Kind alphabet
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE NE_Kind_IntegerAlphabet(MyInteger,MyString)
IMPLICIT NONE
! Passed in variables
DIF3D_Int  MyInteger    ! Some integer
DIF3D_Char MyString     ! The string to return where position (1:1) is always set as "K"
! Local
DIF3D_Int I,K
DIF3D_Int CopyInteger
DIF3D_Int idivisor

idivisor = NE_Kind_SizeAlphabet*NE_Kind_SizeAlphabet*NE_Kind_SizeAlphabet*NE_Kind_SizeAlphabet*NE_Kind_SizeAlphabet
IF ((MyInteger .LT. 0) .OR. (MyInteger .GT. idivisor))  THEN
   WRITE(MODULE_OUT,'("[NE_Kind]...Invalid integer size in IntegerAlphabet?!? ",I9)') MyInteger
   CALL Basic_Abort
END IF

MyString = ""
MyString(1:1) = "K"
CopyInteger = MyInteger
DO I = 1,5
   idivisor = idivisor / NE_Kind_SizeAlphabet
   K = CopyInteger/idivisor
   CopyInteger = CopyInteger - K*idivisor
   MyString(I+1:I+1) = NE_Kind_MyAlphabet(K+1:K+1)
END DO
MyString = ADJUSTL(MyString)                ! Left justify
MyString = TRIM(MyString)                   ! Trim

END SUBROUTINE NE_Kind_IntegerAlphabet

!---------------------------------------------------------------------------------------------------------------------------------
!  NE_Kind_AlphabetInteger(MyInteger,MyString)
!  Takes an 6 letter NE_Kind alphabet word and converts it to an integer using 
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE NE_Kind_AlphabetInteger(MyInteger,MyString)
IMPLICIT NONE
! Passed in variables
DIF3D_Int  MyInteger    ! The integer to return
DIF3D_Char MyString     ! The string to convert where position (1:1) is always ignored
! Local
DIF3D_Int I,J,K
DIF3D_Int idivisor

MyInteger = 0
idivisor = 1
DO I = 5,1,-1
   J = 0
   DO K = 1,NE_Kind_SizeAlphabet
      IF (NE_Kind_MyAlphabet(K:K) .EQ. MyString(I+1:I+1)) THEN
         J = K ! The matching alphabet letter
         EXIT
      END IF
   END DO
   IF (J .EQ. 0) THEN
      WRITE(MODULE_OUT,'("[NE_Kind]...Encountered an invalid character?!? ",A1)') MyString(I+1:I+1)
      CALL Basic_Abort
   END IF
   MyInteger = MyInteger + (J-1)*idivisor
   idivisor = idivisor * NE_Kind_SizeAlphabet
END DO

END SUBROUTINE NE_Kind_AlphabetInteger

END MODULE NE_Kind
