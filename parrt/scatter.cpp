#include "scatter.h"
# include <stdlib.h>
# include <stdio.h>
# include <iostream>
# include <iomanip>

//comm to cells
void comm2cells(int ndof_per_cell, int ncells_per_comm, PetscScalar *u, PetscScalar **v, int forward)
{

  for(int i=0;i<ncells_per_comm;i++)
    for(int j=0;j<ndof_per_cell;j++)
      if(forward>0)
	v[i][j]=u[j+i*ndof_per_cell];
      else
	u[j+i*ndof_per_cell]=v[i][j];
}

void comm2cells(int ndof_per_cell, int ncells_per_comm, PetscReal *u, PetscReal **v, int forward)
{

  for(int i=0;i<ncells_per_comm;i++)
    for(int j=0;j<ndof_per_cell;j++)
      if(forward)
	v[i][j]=u[j+i*ndof_per_cell];
      else
	u[j+i*ndof_per_cell]=v[i][j];

}

//tandem to comm
void tandem2comm(int colour, int ndof_per_comm, PetscReal *dof_global, PetscReal *dof_local)
{

  for(int i=0;i<ndof_per_comm;i++)
    dof_local[i] = dof_global[i+colour*ndof_per_comm];
  
}

//comm to tandem
void comm2tandem(MPI_Comm comm_global, MPI_Comm comm_local, int colour, int ndof_per_comm, int ncomms, PetscReal *dof_local, PetscReal *dof_global)
{
  int ndof=ndof_per_comm * ncomms;
  PetscReal *tmp = (PetscReal *) malloc(ndof*sizeof(PetscReal));

  int rank_local;
  MPI_Comm_rank(comm_local,&rank_local);
  for(int i=0;i<ndof;i++){
    int i_local = i-colour*ndof_per_comm;
    if(rank_local==0 && 0 <= i_local && i_local < ndof_per_comm)
      tmp[i]=dof_local[i_local];
    else
      tmp[i]=0.0;
  }

  MPI_Allreduce(tmp,dof_global,ndof,MPI_DOUBLE,MPI_SUM,comm_global);

  free(tmp);

}

//tandem to comm
void tandem2comm(int colour, int ndof_per_comm, PetscScalar *dof_global, PetscScalar *dof_local)
{

  for(int i=0;i<ndof_per_comm;i++)
    dof_local[i] = dof_global[i+colour*ndof_per_comm];
  
}

//comm to tandem
void comm2tandem(MPI_Comm comm_global, MPI_Comm comm_local, int colour, int ndof_per_comm, int ncomms, PetscScalar *dof_local, PetscScalar *dof_global)
{
  int ndof=ndof_per_comm * ncomms;
  PetscScalar *tmp = (PetscScalar *) malloc(ndof*sizeof(PetscScalar));

  int rank_local;
  MPI_Comm_rank(comm_local,&rank_local);
  for(int i=0;i<ndof;i++){
    int i_local = i-colour*ndof_per_comm;
    if(rank_local==0 && 0 <= i_local && i_local < ndof_per_comm)
      tmp[i]=dof_local[i_local];
    else
      tmp[i]=0.0 + PETSC_i*0.0;
  }

  MPI_Allreduce(tmp,dof_global,ndof,MPIU_SCALAR,MPI_SUM,comm_global);

  free(tmp);

}

void tandem2natural(PetscReal *u, PetscReal *v, int nx, int ny, int nlayers, int ncells_x, int ncells_y, int forward)
{

  int Nx=nx*ncells_x;
  int Ny=ny*ncells_y;

  for(int ic_y=0;ic_y<ncells_y;ic_y++){
    for(int ic_x=0;ic_x<ncells_x;ic_x++){
      for(int il=0;il<nlayers;il++){
	for(int iy=0;iy<ny;iy++){
	  for(int ix=0;ix<nx;ix++){

	    int i=ix+nx*iy+nx*ny*il+nx*ny*nlayers*ic_x+nx*ny*nlayers*ncells_x*ic_y;
	    
	    int jx=ix+ic_x*nx;
	    int jy=iy+ic_y*ny;
	    int jl=il;
	    int j = jx + Nx*jy + Nx*Ny*jl;

	    if(forward>0)
	      v[j]=u[i];
	    else
	      u[i]=v[j];
	      
	  }
	}
      }
    }
  }

}

void tandem2natural(PetscScalar *u, PetscScalar *v, int nx, int ny, int nlayers, int ncells_x, int ncells_y, int forward)
{

  int Nx=nx*ncells_x;
  int Ny=ny*ncells_y;

  for(int ic_y=0;ic_y<ncells_y;ic_y++){
    for(int ic_x=0;ic_x<ncells_x;ic_x++){
      for(int il=0;il<nlayers;il++){
	for(int iy=0;iy<ny;iy++){
	  for(int ix=0;ix<nx;ix++){

	    int i=ix+nx*iy+nx*ny*il+nx*ny*nlayers*ic_x+nx*ny*nlayers*ncells_x*ic_y;
	    
	    int jx=ix+ic_x*nx;
	    int jy=iy+ic_y*ny;
	    int jl=il;
	    int j = jx + Nx*jy + Nx*Ny*jl;

	    if(forward>0)
	      v[j]=u[i];
	    else
	      u[i]=v[j];
	      
	  }
	}
      }
    }
  }

}
