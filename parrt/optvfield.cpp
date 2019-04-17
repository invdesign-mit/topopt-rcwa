# include "optvfield.h"
# include <stdlib.h>
# include <stdio.h>
# include <iostream>
# include <iomanip>

extern int count;

double vdoteglobal(int ndof, double *dof, double *grad, void *data)
{

  lstcell_ * ptdata=(lstcell_ *) data;

  MPI_Comm comm=ptdata->comm;
  int colour=ptdata->colour;
  int ndof_per_cell=ptdata->ndof_per_cell;
  int ncells_per_comm=ptdata->ncells_per_comm;
  int ncomms=ptdata->ncomms;
  cell_ **cell=ptdata->cell;
  FiltersToolBox *flt=ptdata->flt;

  int nlayers=ptdata->nlayers;
  int nlayers_patterned=ptdata->nlayers_patterned;
  double total_thickness=ptdata->total_thickness;
  double *t_min=ptdata->t_min;
  double *t_max=ptdata->t_max;
    
  PetscScalar *dfxdt_per_cell=(PetscScalar *) malloc(2*nlayers_patterned*sizeof(PetscScalar));
  PetscScalar *dfydt_per_cell=(PetscScalar *) malloc(2*nlayers_patterned*sizeof(PetscScalar));
  PetscScalar *dfzdt_per_cell=(PetscScalar *) malloc(2*nlayers_patterned*sizeof(PetscScalar));
  PetscScalar *dfxdt_per_comm=(PetscScalar *) malloc((nlayers-1)*sizeof(PetscScalar));
  PetscScalar *dfydt_per_comm=(PetscScalar *) malloc((nlayers-1)*sizeof(PetscScalar));
  PetscScalar *dfzdt_per_comm=(PetscScalar *) malloc((nlayers-1)*sizeof(PetscScalar));
  PetscScalar *dfxdt_global=(PetscScalar *) malloc((nlayers-1)*sizeof(PetscScalar));
  PetscScalar *dfydt_global=(PetscScalar *) malloc((nlayers-1)*sizeof(PetscScalar));
  PetscScalar *dfzdt_global=(PetscScalar *) malloc((nlayers-1)*sizeof(PetscScalar));
  if(2*nlayers_patterned!=nlayers-1) SETERRQ(comm,1,"2*nlayers_patterned != nlayers-1");

  for(int i=0;i<2*nlayers_patterned;i++){
    dfxdt_per_comm[i]=0.0;
    dfydt_per_comm[i]=0.0;
    dfzdt_per_comm[i]=0.0;
  }
  int ndof_eps=ndof-2*nlayers_patterned;
  
  int ndof_per_comm=ndof_per_cell*ncells_per_comm;
  if(ndof_per_comm*ncomms!=ndof_eps) SETERRQ(comm,1,"ndof_eps != ndof_per_comm * ncomms");

  PetscScalar vxdote_cplx_per_comm=0, vydote_cplx_per_comm=0, vzdote_cplx_per_comm=0;
  PetscReal *filtered_dofreal_per_cell=(PetscReal *) malloc(ndof_per_cell*sizeof(PetscReal));
  PetscScalar *filtered_vxgrad_cplx_per_comm=(PetscScalar *)malloc(ndof_per_comm*sizeof(PetscScalar));
  PetscScalar *filtered_vygrad_cplx_per_comm=(PetscScalar *)malloc(ndof_per_comm*sizeof(PetscScalar));
  PetscScalar *filtered_vzgrad_cplx_per_comm=(PetscScalar *)malloc(ndof_per_comm*sizeof(PetscScalar));
  for(int i=0;i<ncells_per_comm;i++){

    double tmp_t=0;
    for(int j=0;j<nlayers-1;j++){
      cell[i]->thickness[j]=t_min[j]+dof[j+ndof_eps]*(t_max[j]-t_min[j]);
      tmp_t+=cell[i]->thickness[j];
    }
    cell[i]->thickness[nlayers-1]=total_thickness-tmp_t;
    
    filters_apply(comm,&dof[colour*ndof_per_comm + i*ndof_per_cell],filtered_dofreal_per_cell,flt,1);

    double uobj[2],vobj[2],wobj[2];
    v_dot_E(filtered_dofreal_per_cell, cell[i],
	    cell[i]->bx, cell[i]->by,
	    cell[i]->vx_cx, cell[i]->vx_cy, uobj, cell[i]->vxgrad_complex, dfxdt_per_cell,
	    cell[i]->vy_cx, cell[i]->vy_cy, vobj, cell[i]->vygrad_complex, dfydt_per_cell,
	    cell[i]->vz_cx, cell[i]->vz_cy, wobj, cell[i]->vzgrad_complex, dfzdt_per_cell);

    vxdote_cplx_per_comm += uobj[0] + PETSC_i * uobj[1];
    vydote_cplx_per_comm += vobj[0] + PETSC_i * vobj[1];
    vzdote_cplx_per_comm += wobj[0] + PETSC_i * wobj[1];
    
    filters_apply(comm,cell[i]->vxgrad_complex,&filtered_vxgrad_cplx_per_comm[i*ndof_per_cell],flt,-1);
    filters_apply(comm,cell[i]->vygrad_complex,&filtered_vygrad_cplx_per_comm[i*ndof_per_cell],flt,-1);
    filters_apply(comm,cell[i]->vzgrad_complex,&filtered_vzgrad_cplx_per_comm[i*ndof_per_cell],flt,-1);

    for(int j1=0;j1<2*nlayers_patterned;j1++){
      PetscScalar tmp_dfxdt = 0.0;
      PetscScalar tmp_dfydt = 0.0;
      PetscScalar tmp_dfzdt = 0.0;
      for(int j2=j1;j2<2*nlayers_patterned;j2++){
	tmp_dfxdt += dfxdt_per_cell[j2];
	tmp_dfydt += dfydt_per_cell[j2];
	tmp_dfzdt += dfzdt_per_cell[j2];
      }
      dfxdt_per_comm[j1] += tmp_dfxdt;
      dfydt_per_comm[j1] += tmp_dfydt;
      dfzdt_per_comm[j1] += tmp_dfzdt;
    }
    
  }
				       
  MPI_Barrier(comm);
  MPI_Barrier(PETSC_COMM_WORLD);

  int rank;
  MPI_Comm_rank(comm, &rank);
  if(rank!=0){
    vxdote_cplx_per_comm=0.0;
    vydote_cplx_per_comm=0.0;
    vzdote_cplx_per_comm=0.0;
  }
  PetscScalar vdote_cplx_per_comm[3]={vxdote_cplx_per_comm,vydote_cplx_per_comm,vzdote_cplx_per_comm};
  PetscScalar vdote_cplx_global[3];
  MPI_Allreduce(vdote_cplx_per_comm,vdote_cplx_global,3,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);
  double Esq_global=std::pow(std::abs(vdote_cplx_global[0]),2) + std::pow(std::abs(vdote_cplx_global[1]),2) + std::pow(std::abs(vdote_cplx_global[2]),2);
  if(ptdata->verbose_global) PetscPrintf(PETSC_COMM_WORLD,"------at step %d, \" sum_xyz | sum_all v_[x,y,z] . E dr |^2 \" [global] = %.16g \n",
					 count,
					 Esq_global);

  double *tmp_grad=(double*) malloc(ndof*sizeof(double));
  for(int i=0;i<ndof;i++){
    int i_local=i-colour*ndof_per_comm;
    if(rank==0 && 0<=i_local && i_local<ndof_per_comm)
      tmp_grad[i] = \
	  2.0*std::real(std::conj(vdote_cplx_global[0])*filtered_vxgrad_cplx_per_comm[i_local]) \
	+ 2.0*std::real(std::conj(vdote_cplx_global[1])*filtered_vygrad_cplx_per_comm[i_local]) \
	+ 2.0*std::real(std::conj(vdote_cplx_global[2])*filtered_vzgrad_cplx_per_comm[i_local]);
    else
      tmp_grad[i] = 0.0;
  }
  MPI_Barrier(comm);
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Allreduce(tmp_grad,grad,ndof,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Barrier(comm);
  MPI_Barrier(PETSC_COMM_WORLD);
  free(tmp_grad);

  if(rank!=0){
    for(int i=0;i<2*nlayers_patterned;i++){
      dfxdt_per_comm[i]=0.0;
      dfydt_per_comm[i]=0.0;
      dfzdt_per_comm[i]=0.0;
    }
  }
  MPI_Allreduce(dfxdt_per_comm,dfxdt_global,nlayers-1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(dfydt_per_comm,dfydt_global,nlayers-1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(dfzdt_per_comm,dfzdt_global,nlayers-1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);
  for(int i=0;i<nlayers-1;i++){
        grad[i+ndof_eps] = \
	  2.0*std::real(std::conj(vdote_cplx_global[0])*dfxdt_global[i]) + \
	  2.0*std::real(std::conj(vdote_cplx_global[1])*dfydt_global[i]) + \
	  2.0*std::real(std::conj(vdote_cplx_global[2])*dfzdt_global[i]);
	grad[i+ndof_eps]=(t_max[i]-t_min[i])*grad[i+ndof_eps];
  }
  
  free(filtered_dofreal_per_cell);
  free(filtered_vxgrad_cplx_per_comm);
  free(filtered_vygrad_cplx_per_comm);
  free(filtered_vzgrad_cplx_per_comm);

  free(dfxdt_per_cell);
  free(dfxdt_per_comm);
  free(dfxdt_global);
  free(dfydt_per_cell);
  free(dfydt_per_comm);
  free(dfydt_global);
  free(dfzdt_per_cell);
  free(dfzdt_per_comm);
  free(dfzdt_global);
  
  int print_at=ptdata->print_at;
  if(count%print_at==0){

    int *ngrid=ptdata->ngrid;
    int nlayers_patterned=ptdata->nlayers_patterned;
    int ncells_along_x=ptdata->ncells_along_x;
    int ncells_along_y=ptdata->ncells_along_y;
    int ndof_global=ndof_per_comm*ncomms;

    double *dof_raw_natural=(double *)malloc(ndof_global*sizeof(double));
    double *dof_filtered_comm=(double*)malloc(ndof_per_comm*sizeof(double));
    double *dof_filtered_tandem=(double *)malloc(ndof_global*sizeof(double));
    double *dof_filtered_natural=(double *)malloc(ndof_global*sizeof(double));
    tandem2natural(dof, dof_raw_natural, ngrid[0], ngrid[1], nlayers_patterned, ncells_along_x, ncells_along_y, 1);
    for(int i=0;i<ncells_per_comm;i++)
      filters_apply(comm,&dof[colour*ndof_per_comm+i*ndof_per_cell],&dof_filtered_comm[i*ndof_per_cell],flt,1);
    comm2tandem(PETSC_COMM_WORLD, comm, colour, ndof_per_comm, ncomms, dof_filtered_comm, dof_filtered_tandem);
    tandem2natural(dof_filtered_tandem, dof_filtered_natural, ngrid[0], ngrid[1], nlayers_patterned, ncells_along_x, ncells_along_y, 1);
        
    MPI_Barrier(PETSC_COMM_WORLD);
        
    std::string fname1="step"+std::to_string(count)+"_dof_raw_tandem.dat";
    std::string fname2="step"+std::to_string(count)+"_dof_raw_natural.dat";
    std::string fname3="step"+std::to_string(count)+"_dof_filtered_tandem.dat";
    std::string fname4="step"+std::to_string(count)+"_dof_filtered_natural.dat";
    if(ptdata->print_dof_raw_tandem) writetofiledouble(PETSC_COMM_WORLD,(char *)fname1.c_str(),dof,ndof_global);
    if(ptdata->print_dof_raw_natural) writetofiledouble(PETSC_COMM_WORLD,(char *)fname2.c_str(),dof_raw_natural,ndof_global);
    if(ptdata->print_dof_filtered_tandem) writetofiledouble(PETSC_COMM_WORLD,(char *)fname3.c_str(),dof_filtered_tandem,ndof_global);
    if(ptdata->print_dof_filtered_natural) writetofiledouble(PETSC_COMM_WORLD,(char *)fname4.c_str(),dof_filtered_natural,ndof_global);

    std::string fname5="step"+std::to_string(count)+"_thicknesses.dat";
    if(ptdata->print_dof_thickness) writetofiledouble(PETSC_COMM_WORLD,(char *)fname5.c_str(),&(dof[ndof_global]),nlayers-1);
    
    MPI_Barrier(PETSC_COMM_WORLD);

    free(dof_raw_natural);
    free(dof_filtered_comm);
    free(dof_filtered_tandem);
    free(dof_filtered_natural);

  }

  count++;
  return Esq_global;
  
}

