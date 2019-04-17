#ifndef GUARD_output_h
#define GUARD_output_h

#include "petsc.h"

void writetofile(MPI_Comm comm, char *name, PetscScalar *data, PetscInt n);

void writetofiledouble(MPI_Comm comm, char *name, double *data, PetscInt n);

PetscErrorCode saveVecMfile(MPI_Comm comm, Vec vec, const char *filename, const char *dsetname);

#endif
