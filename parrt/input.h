#ifndef GUARD_input_h
#define GUARD_input_h

#include "petsc.h"

#define PI 3.14159265358979323846

PetscErrorCode getreal(const char *flag, double *var, double autoval);

PetscErrorCode getint(const char *flag, int *var, int autoval);

PetscErrorCode getstr(const char *flag, char *filename, const char default_filename[]);

PetscErrorCode getintarray(const char *flag, int *z, int *nz, int default_val);

PetscErrorCode getrealarray(const char *flag, PetscReal *r, int *nr, PetscReal default_val);

void readfromfile(char *name, PetscScalar *data, PetscInt n);

void readfromfiledouble(char *name, double *data, PetscInt n);

#endif
