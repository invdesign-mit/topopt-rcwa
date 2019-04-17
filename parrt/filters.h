#ifndef GUARD_filters_h
#define GUARD_filters_h

# include "petsc.h"
# include "petscsys.h"
# include "output.h"

typedef struct{

  int filter_choice;
  int ndof;
  double filter_threshold;
  double filter_steepness;
  int filter_dim;
  Mat W;
  Vec rho_grad;

}FiltersToolBox;

typedef void (*FilterType)(Vec rho_in, Vec rho_out, FiltersToolBox *flt);

void create_density_filter_1d(MPI_Comm comm, Mat *Wout, int radius_pixels, int ndof_per_layer, int num_layers);

void create_density_filter_2d(MPI_Comm comm, Mat *Wout, PetscReal rp, int mx_per_layer, int my_per_layer, int num_layers);

void threshold_projection_filter(Vec rho_in, Vec rho_out, Vec rho_grad, double filter_threshold, double filter_steepness);

PetscErrorCode array2mpi(PetscScalar *pt, Vec v);
PetscErrorCode array2mpi(PetscReal *pt, Vec v);

PetscErrorCode mpi2array(Vec v, PetscScalar *pt, int n);
PetscErrorCode mpi2array(Vec v, PetscReal *pt, int n);

void smear_projection(Vec rho_in, Vec rho_out, FiltersToolBox *flt);

void smear_projection_undo(Vec rho_in, Vec rho_out, FiltersToolBox *flt);

void filters_apply(MPI_Comm comm, PetscScalar *u_in, PetscScalar *u_out, FiltersToolBox *flt, int filter_direction);
void filters_apply(MPI_Comm comm, PetscReal *u_in, PetscReal *u_out, FiltersToolBox *flt, int filter_direction);

void filters_initialize(MPI_Comm comm, const int filter_choice, const double filter_threshold, const double filter_steepness, const int filter_dim, const double radius, const int *dofinfo, FiltersToolBox *flt);

#endif
