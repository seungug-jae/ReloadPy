!---------------------------------------------------------------------------------------------------------------------------------
! Assembles the table of timing data from all of the processors and prints it to the screen
!  Copyright(c) 2005 Argonne National Laboratory
!---------------------------------------------------------------------------------------------------------------------------------
#undef __FUNCT__
#define __FUNCT__ "Basic_Abort"
SUBROUTINE Basic_Abort()
IMPLICIT NONE
! MPI preprocessing information handled via PETSc interface
#include "DIF3D_Types.h"

DIF3D_Int ReturnedError,I(10)

! Make certain everyone is ready to quit at the same time
CALL  ABORT

END SUBROUTINE Basic_Abort
