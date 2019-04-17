#ifndef GUARD_multiobj_h
#define GUARD_multiobj_h

#include "petsc.h"
#include "nlopt.h"
#include "filters.h"
#include "output.h"
#include "scatter.h"
#include <string>
#include <sstream>
#include "optvfield.h"

typedef struct{
  
  int *ngrid;
  int nlayers_patterned;
  int ncells_along_x;
  int ncells_along_y;
  int ndof_global;
  int ncomms;
  int ncells_per_comm;
  int ndof_per_cell;
  int colour;
  MPI_Comm comm;
  int print_at;
  int print_dof_raw_tandem;
  int print_dof_raw_natural;
  int print_dof_filtered_tandem;
  int print_dof_filtered_natural;
  int print_dof_thickness;
  FiltersToolBox *flt;
  
} PrintInfo_;

typedef struct{

  int nspecs;
  lstcell_ **lstcell;

} avglist_;

double dummyobj(int ndofAll, double *dofAll, double *gradAll, void *data);

double maximin_vdote(int ndofAll, double *dofAll, double *gradAll, void *data);

#endif
