#ifndef GUARD_scatter_h
#define GUARD_scatter_h

#include "petsc.h"

void comm2cells(int ndof_per_cell, int ncells_per_comm, PetscScalar *u, PetscScalar **v, int forward);
void comm2cells(int ndof_per_cell, int ncells_per_comm, PetscReal *u, PetscReal **v, int forward);

void tandem2comm(int colour, int ndof_per_comm, PetscReal *dof_global, PetscReal *dof_local);
void tandem2comm(int colour, int ndof_per_comm, PetscScalar *dof_global, PetscScalar *dof_local);

void comm2tandem(MPI_Comm comm_global, MPI_Comm comm_local, int colour, int ndof_per_comm, int ncomms, PetscReal *dof_local, PetscReal *dof_global);
void comm2tandem(MPI_Comm comm_global, MPI_Comm comm_local, int colour, int ndof_per_comm, int ncomms, PetscScalar *dof_local, PetscScalar *dof_global);

void tandem2natural(PetscReal *u, PetscReal *v, int nx, int ny, int nlayers, int ncells_x, int ncells_y, int forward);
void tandem2natural(PetscScalar *u, PetscScalar *v, int nx, int ny, int nlayers, int ncells_x, int ncells_y, int forward);

#endif
