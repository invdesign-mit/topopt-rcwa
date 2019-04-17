# include "multiobj.h"
# include "optvfield.h"
# include <stdlib.h>
# include <stdio.h>
# include <iostream>
# include <iomanip>

extern int count;

double dummyobj(int ndofAll, double *dofAll, double *gradAll, void *data)
{

  int ndof=ndofAll-1;
  for(int i=0;i<ndof;i++)
    gradAll[i]=0.0;
  gradAll[ndof]=1.0;
  count++;
  PetscPrintf(PETSC_COMM_WORLD,"---at step %d, dummyobj is: %0.8g\n",count,dofAll[ndof]);
  
  PrintInfo_ *ptdata = (PrintInfo_ *)data;
  int print_at=ptdata->print_at;
  if(count%print_at==0){

    int ndof_per_cell=ptdata->ndof_per_cell;
    int ncells_per_comm=ptdata->ncells_per_comm;
    int ndof_per_comm=ndof_per_cell*ncells_per_comm;
    int colour=ptdata->colour;
    int ncomms=ptdata->ncomms;
    int *ngrid=ptdata->ngrid;
    int nlayers_patterned=ptdata->nlayers_patterned;
    int ncells_along_x=ptdata->ncells_along_x;
    int ncells_along_y=ptdata->ncells_along_y;
    int ndof_global=ndof_per_comm*ncomms;
    FiltersToolBox *flt=ptdata->flt;
    MPI_Comm comm=ptdata->comm;

    double *dof_raw_natural=(double *)malloc(ndof_global*sizeof(double));
    double *dof_filtered_comm=(double*)malloc(ndof_per_comm*sizeof(double));
    double *dof_filtered_tandem=(double *)malloc(ndof_global*sizeof(double));
    double *dof_filtered_natural=(double *)malloc(ndof_global*sizeof(double));
    tandem2natural(&(dofAll[0]), dof_raw_natural, ngrid[0], ngrid[1], nlayers_patterned, ncells_along_x, ncells_along_y, 1);
    for(int i=0;i<ncells_per_comm;i++)
      filters_apply(comm,&(dofAll[colour*ndof_per_comm+i*ndof_per_cell]),&dof_filtered_comm[i*ndof_per_cell],flt,1);
    comm2tandem(PETSC_COMM_WORLD, comm, colour, ndof_per_comm, ncomms, dof_filtered_comm, dof_filtered_tandem);
    tandem2natural(dof_filtered_tandem, dof_filtered_natural, ngrid[0], ngrid[1], nlayers_patterned, ncells_along_x, ncells_along_y, 1);

    MPI_Barrier(PETSC_COMM_WORLD);

    std::string fname1="step"+std::to_string(count)+"_dof_raw_tandem.dat";
    std::string fname2="step"+std::to_string(count)+"_dof_raw_natural.dat";
    std::string fname3="step"+std::to_string(count)+"_dof_filtered_tandem.dat";
    std::string fname4="step"+std::to_string(count)+"_dof_filtered_natural.dat";
    if(ptdata->print_dof_raw_tandem) writetofiledouble(PETSC_COMM_WORLD,(char *)fname1.c_str(),&(dofAll[0]),ndof_global);
    if(ptdata->print_dof_raw_natural) writetofiledouble(PETSC_COMM_WORLD,(char *)fname2.c_str(),dof_raw_natural,ndof_global);
    if(ptdata->print_dof_filtered_tandem) writetofiledouble(PETSC_COMM_WORLD,(char *)fname3.c_str(),dof_filtered_tandem,ndof_global);
    if(ptdata->print_dof_filtered_natural) writetofiledouble(PETSC_COMM_WORLD,(char *)fname4.c_str(),dof_filtered_natural,ndof_global);

    std::string fname5="step"+std::to_string(count)+"_thicknesses.dat";
    if(ptdata->print_dof_thickness) writetofiledouble(PETSC_COMM_WORLD,(char *)fname5.c_str(),&(dofAll[ndof_global]),ndofAll-ndof_global-1);
    
    MPI_Barrier(PETSC_COMM_WORLD);

    free(dof_raw_natural);
    free(dof_filtered_comm);
    free(dof_filtered_tandem);
    free(dof_filtered_natural);

  }
  
  return dofAll[ndof];
  
}


double maximin_vdote(int ndofAll, double *dofAll, double *gradAll, void *data)
{

  int ndof=ndofAll-1;
  double objval=vdoteglobal(ndof,&(dofAll[0]),&(gradAll[0]),data);
  count--;
  for(int i=0;i<ndof;i++)
    gradAll[i]=-gradAll[i];
  gradAll[ndof]=1.0;
  MPI_Barrier(PETSC_COMM_WORLD);
  lstcell_ *ptdata=(lstcell_ *)data;
  PetscPrintf(PETSC_COMM_WORLD,"---at step %d, specID %d, vdoEobj is: %0.16g \n",count,ptdata->specID,objval);
  
  return dofAll[ndof]-objval;

}


