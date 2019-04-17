#ifndef GUARD_optvfield_h
#define GUARD_optvfield_h

#include "petsc.h"
#include "nlopt.h"
#include "objfuncs.h"
#include "wrap.h"
#include "filters.h"
#include "output.h"
#include "scatter.h"
#include <string>
#include <sstream>

typedef struct{

  MPI_Comm comm;
  int colour;
  int ndof_per_cell;
  int ncells_per_comm;
  int ncomms;
  cell_ **cell;
  FiltersToolBox *flt;
  int verbose_global;
  int verbose_comm;
  int verbose_cell;

  int specID;
  
  int print_at;
  int print_dof_raw_tandem;
  int print_dof_raw_natural;
  int print_dof_filtered_tandem;
  int print_dof_filtered_natural;
  int ngrid[2];
  int nlayers_patterned;
  int ncells_along_x;
  int ncells_along_y;

  int ndof_global;
  int nspecs;

  int nlayers;
  double total_thickness;
  double *t_min;
  double *t_max;
  int print_dof_thickness;
  
} lstcell_;

double vdoteglobal(int ndof, double *dof, double *grad, void *data);

#endif
