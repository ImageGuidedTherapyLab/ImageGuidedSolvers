static char help[] = "create dense matrix and multiply \n\n";
/* 
  Include "petscmat.h" so that we can use matrices.
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h    - vectors
     petscmat.h    - matrices
     petscis.h     - index sets            petscviewer.h - viewers               
*/
// C include 
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
// C++ include 
#include <iostream>
// petsc include
#include "petscmat.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  PetscErrorCode ierr;
  Mat                   A,B,C;                /* matrix */

  PetscInitialize(&argc,&args,(char *)0,help);

  ierr=PetscPrintf(PETSC_COMM_WORLD, "export MPIRUN_RANK=%s  \n",getenv("MPIRUN_RANK"  ) );CHKERRQ(ierr);
  ierr=PetscPrintf(PETSC_COMM_WORLD, "export MPIRUN_NPROCS=%s\n",getenv("MPIRUN_NPROCS") );CHKERRQ(ierr);
  ierr=PetscPrintf(PETSC_COMM_WORLD, "export MPIRUN_ID=%s    \n",getenv("MPIRUN_ID"    ) );CHKERRQ(ierr);
  ierr=PetscPrintf(PETSC_COMM_WORLD, "export MPIRUN_HOST=%s  \n",getenv("MPIRUN_HOST"  ) );CHKERRQ(ierr);
  ierr=PetscPrintf(PETSC_COMM_WORLD, "export MPIRUN_PORT=%s  \n",getenv("MPIRUN_PORT"  ) );CHKERRQ(ierr);

  PetscInt nsize = 1000;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-nsize",&nsize,PETSC_NULL);
  ierr=PetscPrintf(PETSC_COMM_WORLD, "nsize %d...\n",nsize);CHKERRQ(ierr);
  ierr = MatCreateMPIDense(PETSC_COMM_WORLD,nsize,nsize,nsize,nsize,
                                   PETSC_NULL,&A);CHKERRQ(ierr);
  ierr = MatCreateMPIDense(PETSC_COMM_WORLD,nsize,nsize,nsize,nsize,
                                   PETSC_NULL,&B);CHKERRQ(ierr);

  ierr=PetscPrintf(PETSC_COMM_WORLD, "multiply...\n");CHKERRQ(ierr);

  ierr = MatMatMult(A,B,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&C); CHKERRQ(ierr);

  ierr=PetscPrintf(PETSC_COMM_WORLD, "cleanup...\n");CHKERRQ(ierr);
  ierr = MatDestroy(A);CHKERRQ(ierr);
  ierr = MatDestroy(B);CHKERRQ(ierr);
  ierr = MatDestroy(C);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

