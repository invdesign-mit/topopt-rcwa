# include "rcwa.h"
# include <stdlib.h>
# include <stdio.h>
# include <limits>
# include <math.h>
# include <complex>
# include "gsel.h"
# include "fft_iface.h"
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <string>
# include <sstream>
# include "numalloc.h"
# include <TBLAS.h>
#ifdef HAVE_BLAS
# include <TBLAS_ext.h>
#endif
# include "mpi.h"
# include "petsc.h"
# include "petscsys.h"
# include "filters.h"
# include "input.h"
# include "output.h"
# include "wrap.h"
# include "objfuncs.h"
# include "optitemp.h"
# include "optvfield.h"
# include "scatter.h"
# include "near2far.h"
# include "specs.h"
# include "multiobj.h"

int count=0;

int main(int argc, char *argv[]){

  MPI_Init(NULL, NULL);
  PetscInitialize(&argc,&argv,NULL,NULL);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  PetscPrintf(PETSC_COMM_WORLD,"\tThe total number of processors is %d\n",size);

  int ncells_along_x, ncells_along_y, ncells_per_comm;
  getint("-ncells_along_x",&ncells_along_x,1);
  getint("-ncells_along_y",&ncells_along_y,1);
  getint("-ncells_per_comm",&ncells_per_comm,1);

  int ncells,ncomms;
  ncells=ncells_along_x*ncells_along_y;
  ncomms=ncells/ncells_per_comm;
  if(!(size%ncomms==0)) SETERRQ(PETSC_COMM_WORLD,1,"The number of processes must be divisible by the number of subcomms.");

  int np_per_comm=size/ncomms;

  PetscPrintf(PETSC_COMM_WORLD,
	      "\tThe total number of cells is %d: %d cells along x and %d cells along y.\n\tThe number of subcomms is %d, each with %d core.\n\tEach subcomm sequentially handles %d cells.\n",
	      ncells,
	      ncells_along_x,
	      ncells_along_y,
	      ncomms,
	      np_per_comm,
	      ncells_per_comm);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm subcomm;
  int colour = rank/np_per_comm;
  MPI_Comm_split(MPI_COMM_WORLD, colour, rank, &subcomm);
  
  int ngrid[2];
  getint("-ngrid[0]",&ngrid[0],100);
  getint("-ngrid[1]",&ngrid[1],100);
  
  int nlayers_uniform,nlayers_patterned,nG,verbose_local;
  double Lx,Ly;
  getint("-nlayers_uniform",&nlayers_uniform,2);
  getint("-nlayers_patterned",&nlayers_patterned,1);
  getint("-nG",&nG,100);
  getreal("-Lx",&Lx,1.0);
  getreal("-Ly",&Ly,1.0);
  getint("-verbose_local", &verbose_local, 0);
  int nlayers = nlayers_uniform + nlayers_patterned;
  
  int layerID_uniform[nlayers_uniform], layerID_patterned[nlayers_patterned],nz_integrate_patterned[nlayers_patterned];
  double thickness[nlayers];
  getintarray("-uniform_layers",   layerID_uniform,   &nlayers_uniform,0);
  getintarray("-patterned_layers", layerID_patterned, &nlayers_patterned,1);
  getintarray("-nz",               nz_integrate_patterned,      &nlayers_patterned,100);
  getrealarray("-thickness",       thickness,         &nlayers,0.0);

  double total_thickness;
  double *t_min=(double *)malloc((nlayers-1)*sizeof(double));
  double *t_max=(double *)malloc((nlayers-1)*sizeof(double));
  int print_dof_thickness;
  getreal("-total_thickness", &total_thickness, 4);
  int nget=nlayers-1;
  getrealarray("-t_min", t_min, &nget, 1);
  getrealarray("-t_max", t_max, &nget, 1);
  getint("-print_dof_thickness", &print_dof_thickness, 1);
  
  double oxy[2];
  nget=2;
  getrealarray("-origin_xy",oxy,&nget,0.0);

  
  int Ngrid[2]={ngrid[0]*ncells_along_x,ngrid[1]*ncells_along_y};
  int ng2_per_cell=ngrid[0]*ngrid[1];
  int ng2_per_comm=ng2_per_cell*ncells_per_comm;
  int ng2_global=Ngrid[0]*Ngrid[1];
  int ndof_per_cell=ng2_per_cell*nlayers_patterned;
  int ndof_per_comm=ndof_per_cell*ncells_per_comm;
  int ndof_global=ndof_per_comm*ncomms;

  PetscPrintf(PETSC_COMM_WORLD,"****Attention: the ispec_float must be given in the order freq,polar,azimuth,input_polarization[4],xyz_far[3].\n");
  int nspecs;
  getint("-nspecs",&nspecs,1);
  cell_ *cell = (cell_*)malloc(nspecs*ncells_per_comm*sizeof(cell_));
  for(int ispec=0;ispec<nspecs;ispec++){

    std::string spec_name="-ispec"+std::to_string(ispec)+"_float";
    double specs_float[10];
    nget=10;
    getrealarray((char *)spec_name.c_str(),specs_float,&nget,1.0);
    double freq=specs_float[0];
    double polar=specs_float[1];
    double azimuth=specs_float[2];
    double input_polarization[4]={specs_float[3],specs_float[4],specs_float[5],specs_float[6]};
    double xyz_far[3]={specs_float[7],specs_float[8],specs_float[9]};

    spec_name="-ispec"+std::to_string(ispec)+"_symxy";
    int symxy[2];
    nget=2;
    getintarray((char *)spec_name.c_str(),symxy,&nget,0);

    double n_uniform[nlayers_uniform],epsdiff_patterned[nlayers_patterned],epsbkg_patterned[nlayers_patterned];
    spec_name="-ispec"+std::to_string(ispec)+"_nuniform";
    getrealarray((char *)spec_name.c_str(), n_uniform,         &nlayers_uniform,1.0);
    spec_name="-ispec"+std::to_string(ispec)+"_epsdiff";
    getrealarray((char *)spec_name.c_str(), epsdiff_patterned, &nlayers_patterned,2.0);
    spec_name="-ispec"+std::to_string(ispec)+"_epsbkg";
    getrealarray((char *)spec_name.c_str(), epsbkg_patterned,  &nlayers_patterned,1.0);

    for(int icell=0;icell<ncells_per_comm;icell++){

      int i=icell+ispec*ncells_per_comm;

      cell[i].nlayers_uniform=nlayers_uniform;
      cell[i].nlayers_patterned=nlayers_patterned;
      cell[i].nG=(unsigned int)std::abs(nG);
      cell[i].Lr[0]=Lx;
      cell[i].Lr[3]=Ly;
      cell[i].freq=freq;
      cell[i].dir_cosine[0] = n_uniform[0] * std::sin(polar * M_PI/180) * std::cos(azimuth * M_PI/180);
      cell[i].dir_cosine[1] = n_uniform[0] * std::sin(polar * M_PI/180) * std::sin(azimuth * M_PI/180);
      cell[i].ngrid[0]=ngrid[0];
      cell[i].ngrid[1]=ngrid[1];
      
      set_cell(cell+i,
	       layerID_uniform, n_uniform,
	       layerID_patterned, epsdiff_patterned, epsbkg_patterned, nz_integrate_patterned,
	       thickness);
      
    }

    MPI_Barrier(PETSC_COMM_WORLD);
    populate_b_and_c(colour,ngrid,ng2_global,ng2_per_comm,ncells_per_comm,
		     Lx,Ly,ncells_along_x,ncells_along_y,Ngrid,
		     freq,n_uniform[0],polar,azimuth,input_polarization,
		     oxy,symxy,xyz_far,
		     &(cell[ispec*ncells_per_comm]));
    MPI_Barrier(PETSC_COMM_WORLD);

  }
  
  MPI_Barrier(PETSC_COMM_WORLD);
  
  FiltersToolBox flt;
  int dofinfo[3]={ngrid[0],ngrid[1],nlayers_patterned};
  int filter_choice=0, filter_dim=2;
  double filter_threshold, filter_steepness, filter_radius;
  getreal("-filter_threshold",&filter_threshold,0.5);
  getreal("-filter_steepness",&filter_steepness,0.0);
  getreal("-filter_radius",&filter_radius,1.0);
  filters_initialize(subcomm, filter_choice, filter_threshold, filter_steepness, filter_dim, filter_radius, dofinfo, &flt);

  int verbose_global,verbose_comm,verbose_cell;
  int print_at;
  int print_dof_raw_tandem,print_dof_raw_natural;
  int print_dof_filtered_tandem,print_dof_filtered_natural;
  getint("-verbose_global",&(verbose_global),1);
  getint("-verbose_comm",&(verbose_comm),0);
  getint("-verbose_cell",&(verbose_cell),0);
  getint("-print_at",&(print_at),500);
  getint("-print_dof_raw_tandem",&(print_dof_raw_tandem),1);
  getint("-print_dof_raw_natural",&(print_dof_raw_natural),0);
  getint("-print_dof_filtered_tandem",&(print_dof_filtered_tandem),0);
  getint("-print_dof_filtered_natural",&(print_dof_filtered_natural),1);

  lstcell_* lstcell = (lstcell_ *)malloc(nspecs*sizeof(lstcell_));
  for(int j=0;j<nspecs;j++){
    lstcell[j].specID=j;
    lstcell[j].comm=subcomm;
    lstcell[j].colour=colour;
    lstcell[j].ndof_per_cell=ndof_per_cell;
    lstcell[j].ncells_per_comm=ncells_per_comm;
    lstcell[j].ncomms=ncomms;
    lstcell[j].flt=&flt;
    lstcell[j].cell=(cell_ **)malloc(ncells_per_comm*sizeof(cell_*));

    lstcell[j].verbose_global=verbose_global;
    lstcell[j].verbose_comm=verbose_comm;
    lstcell[j].verbose_cell=verbose_cell;
    lstcell[j].print_at=print_at;
    lstcell[j].print_dof_raw_tandem=print_dof_raw_tandem;
    lstcell[j].print_dof_raw_natural=print_dof_raw_natural;
    lstcell[j].print_dof_filtered_tandem=print_dof_filtered_tandem;
    lstcell[j].print_dof_filtered_natural=print_dof_filtered_natural;

    lstcell[j].ngrid[0]=ngrid[0];
    lstcell[j].ngrid[1]=ngrid[1];
    lstcell[j].nlayers_patterned=nlayers_patterned;
    lstcell[j].ncells_along_x=ncells_along_x;
    lstcell[j].ncells_along_y=ncells_along_y;

    for(int i=0;i<ncells_per_comm;i++)
      lstcell[j].cell[i]=&(cell[i+ncells_per_comm*j]);

    lstcell[j].nlayers=nlayers;
    lstcell[j].total_thickness=total_thickness;
    lstcell[j].t_min=(double *)malloc((nlayers-1)*sizeof(double));
    lstcell[j].t_max=(double *)malloc((nlayers-1)*sizeof(double));
    for(int i=0;i<nlayers-1;i++){
      lstcell[j].t_min[i]=t_min[i];
      lstcell[j].t_max[i]=t_max[i];
    }
    lstcell[j].print_dof_thickness=print_dof_thickness;
    
  }

  if(ndof_global!=ng2_global*nlayers_patterned) PetscPrintf(PETSC_COMM_WORLD,"Error: check the TOTAL number of DOFs\n");

  int ndof_tot = ndof_global + nlayers - 1;
  
  int Job;
  getint("-Job",&Job,1);

  if(Job==0){

    int specID;
    getint("-specID",&specID,0);

    double *dof_tandem =  (double *) malloc(ndof_tot*sizeof(double));
    double *grad_tandem = (double *) malloc(ndof_tot*sizeof(double));
    readfromfiledouble((char *)"dof.txt",dof_tandem,ndof_tot);
    MPI_Barrier(PETSC_COMM_WORLD);
    
    int sp;
    double s0,s1,ds;
    getint("-sp",&sp,ndof_global/2);
    getreal("-s0",&s0,0.0);
    getreal("-s1",&s1,1.0);
    getreal("-ds",&ds,0.01);
    for(double s=s0;s<s1;s+=ds){
      dof_tandem[sp]=s;
      double objval=vdoteglobal(ndof_tot,dof_tandem,grad_tandem,&(lstcell[specID]));
      PetscPrintf(PETSC_COMM_WORLD,"objval: %g, %0.16g, %0.16g\n",colour,dof_tandem[sp],objval,grad_tandem[sp]);
    }

    free(dof_tandem);
    free(grad_tandem);

  }
     
  if(Job==1){

    int ndofAll=ndof_tot+1;
    double *dofAll=(double *)malloc(ndofAll*sizeof(double));
    double *lbAll=(double *)malloc(ndofAll*sizeof(double));
    double *ubAll=(double *)malloc(ndofAll*sizeof(double));
    readfromfiledouble((char *)"dof.txt",&(dofAll[0]),ndof_tot);
    getreal("-initial_dummy",&(dofAll[ndofAll-1]),0.0);
    for(int i=0;i<ndof_tot;i++){
      lbAll[i]=0.0;
      ubAll[i]=1.0;
    }
    lbAll[ndofAll-1]=0.0;
    ubAll[ndofAll-1]=std::numeric_limits<double>::infinity();
    MPI_Barrier(PETSC_COMM_WORLD);

    int algouter, alginner, algmaxeval;
    getint("-algouter",&algouter,41);
    getint("-alginner",&alginner,40);
    getint("-algmaxeval",&algmaxeval,500);
    alg_ alg={(nlopt_algorithm)algouter,(nlopt_algorithm)alginner,algmaxeval,1000000,1};

    nlopt_result nlopt_return;

    PrintInfo_ printinfo;
    printinfo.ngrid=ngrid;
    printinfo.nlayers_patterned=nlayers_patterned;
    printinfo.ncells_along_x=ncells_along_x;
    printinfo.ncells_along_y=ncells_along_y;
    printinfo.ndof_global=ndof_global;
    printinfo.ncomms=ncomms;
    printinfo.ncells_per_comm=ncells_per_comm;
    printinfo.ndof_per_cell=ndof_per_cell;
    printinfo.colour=colour;
    printinfo.comm=subcomm;
    printinfo.print_at=print_at;
    printinfo.print_dof_raw_tandem=print_dof_raw_tandem;
    printinfo.print_dof_raw_natural=print_dof_raw_natural;
    printinfo.print_dof_filtered_tandem=print_dof_filtered_tandem;
    printinfo.print_dof_filtered_natural=print_dof_filtered_natural;
    printinfo.print_dof_thickness=print_dof_thickness;
    printinfo.flt=&flt;

    void *constrdata[nspecs];
    nlopt_func* maximins=(nlopt_func*)malloc(nspecs*sizeof(nlopt_func));
    for(int i=0;i<nspecs;i++){
      maximins[i]=(nlopt_func)maximin_vdote;
      lstcell[i].print_at=10*algmaxeval; //since the intermediate outputing is done from the dummyobj function, we should shut down the output from the maximin functions.
      constrdata[i]=&(lstcell[i]);
    }

    double result=optimize_generic(ndofAll, dofAll,
				   lbAll, ubAll,
				   (nlopt_func)dummyobj, &printinfo,
				   maximins,constrdata,nspecs,
				   alg,
				   &nlopt_return);

    PetscPrintf(PETSC_COMM_WORLD,"nlopt return value: %d \n",nlopt_return);
    PetscPrintf(PETSC_COMM_WORLD,"optimal objval: %0.8g \n",result);

    free(dofAll);
    free(lbAll);
    free(ubAll);

  }

  if(Job==2){

    double *dof_tandem=(double *)malloc(ndof_tot*sizeof(double));
    double *lb=(double *)malloc(ndof_tot*sizeof(double));
    double *ub=(double *)malloc(ndof_tot*sizeof(double));
    readfromfiledouble((char *)"dof.txt",dof_tandem,ndof_tot);
    for(int i=0;i<ndof_tot;i++){
      lb[i]=0.0;
      ub[i]=1.0;
    }
    MPI_Barrier(PETSC_COMM_WORLD);

    int algouter, alginner, algmaxeval;
    getint("-algouter",&algouter,41);
    getint("-alginner",&alginner,40);
    getint("-algmaxeval",&algmaxeval,500);
    alg_ alg={(nlopt_algorithm)algouter,(nlopt_algorithm)alginner,algmaxeval,1000000,1};

    nlopt_result nlopt_return;

    int specID;
    getint("-specID", &specID, 0);

    double result=optimize_generic(ndof_tot, dof_tandem,
				   lb, ub,
				   (nlopt_func)vdoteglobal, &(lstcell[specID]),
				   NULL,NULL,0,
				   alg,
				   &nlopt_return);

    PetscPrintf(PETSC_COMM_WORLD,"nlopt return value: %d \n",nlopt_return);
    PetscPrintf(PETSC_COMM_WORLD,"optimal objval: %0.8g \n",result);

    free(dof_tandem);
    free(lb);
    free(ub);

  }

  for(int i=0;i<nspecs*ncells_per_comm;i++)
    free_cell(&(cell[i]));

  MPI_Barrier(PETSC_COMM_WORLD);
  PetscFinalize();

  return 0;
}


