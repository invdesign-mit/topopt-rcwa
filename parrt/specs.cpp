#include "specs.h"

void create_incident_pw(PetscScalar *bx, PetscScalar *by, double freq, double n_medium, double polar, double azimuth, double j[4], double Lx, double Ly, int Ngrid[2])
{

  double k[2];
  k[0]=2 * M_PI * n_medium * freq * std::sin(polar*M_PI/180.0) * std::cos(azimuth*M_PI/180.0);
  k[1]=2 * M_PI * n_medium * freq * std::sin(polar*M_PI/180.0) * std::sin(azimuth*M_PI/180.0);
  double dg[2]={Lx/Ngrid[0],Ly/Ngrid[1]};

  PetscScalar ax=j[0]*std::exp(PETSC_i*j[1]*M_PI/180.0);
  PetscScalar ay=j[2]*std::exp(PETSC_i*j[3]*M_PI/180.0);

  for(int iy=0;iy<Ngrid[1];iy++){
    for(int ix=0;ix<Ngrid[0];ix++){

      int i=ix+Ngrid[0]*iy;
      double x=(ix+0.5)*dg[0];
      double y=(iy+0.5)*dg[1];
      
      bx[i] = ax * std::exp(-PETSC_i * k[0] * x -PETSC_i * k[1] * y ) * (-2.0*M_PI*freq*PETSC_i);
      by[i] = ay * std::exp(-PETSC_i * k[0] * x -PETSC_i * k[1] * y ) * (-2.0*M_PI*freq*PETSC_i);

    }
  }

}

void populate_b_and_c(int colour, int ngrid[2], int ng2_global, int ng2_per_comm, int ncells_per_comm,
		      double Lx, double Ly, int ncells_along_x, int ncells_along_y, int Ngrid[2],
		      double freq, double n_medium, double polar, double azimuth, double input_polarization[4],
		      double oxy[2], int symxy[2], double xyzfar[3],
		      cell_ *cell)
{

  PetscScalar *bx_natural,*bx_tandem,*bx_comm;
  PetscScalar *by_natural,*by_tandem,*by_comm;
  PetscScalar *vxcx_natural,*vxcx_tandem,*vxcx_comm;
  PetscScalar *vycx_natural,*vycx_tandem,*vycx_comm;
  PetscScalar *vzcx_natural,*vzcx_tandem,*vzcx_comm;
  PetscScalar *vxcy_natural,*vxcy_tandem,*vxcy_comm;
  PetscScalar *vycy_natural,*vycy_tandem,*vycy_comm;
  PetscScalar *vzcy_natural,*vzcy_tandem,*vzcy_comm;

  bx_natural=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  bx_tandem=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  bx_comm=(PetscScalar *)malloc(ng2_per_comm*sizeof(PetscScalar));
  by_natural=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  by_tandem=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  by_comm=(PetscScalar *)malloc(ng2_per_comm*sizeof(PetscScalar));
  vxcx_natural=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  vxcx_tandem=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  vxcx_comm=(PetscScalar *)malloc(ng2_per_comm*sizeof(PetscScalar));
  vycx_natural=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  vycx_tandem=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  vycx_comm=(PetscScalar *)malloc(ng2_per_comm*sizeof(PetscScalar));
  vzcx_natural=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  vzcx_tandem=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  vzcx_comm=(PetscScalar *)malloc(ng2_per_comm*sizeof(PetscScalar));
  vxcy_natural=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  vxcy_tandem=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  vxcy_comm=(PetscScalar *)malloc(ng2_per_comm*sizeof(PetscScalar));
  vycy_natural=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  vycy_tandem=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  vycy_comm=(PetscScalar *)malloc(ng2_per_comm*sizeof(PetscScalar));
  vzcy_natural=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  vzcy_tandem=(PetscScalar *)malloc(ng2_global*sizeof(PetscScalar));
  vzcy_comm=(PetscScalar *)malloc(ng2_per_comm*sizeof(PetscScalar));

  int ng2_per_cell=ngrid[0]*ngrid[1];
  PetscScalar **bx_cell=(PetscScalar **)malloc(ncells_per_comm*sizeof(PetscScalar *));
  PetscScalar **by_cell=(PetscScalar **)malloc(ncells_per_comm*sizeof(PetscScalar *));
  PetscScalar **vxcx_cell=(PetscScalar **)malloc(ncells_per_comm*sizeof(PetscScalar *));
  PetscScalar **vycx_cell=(PetscScalar **)malloc(ncells_per_comm*sizeof(PetscScalar *));
  PetscScalar **vzcx_cell=(PetscScalar **)malloc(ncells_per_comm*sizeof(PetscScalar *));
  PetscScalar **vxcy_cell=(PetscScalar **)malloc(ncells_per_comm*sizeof(PetscScalar *));
  PetscScalar **vycy_cell=(PetscScalar **)malloc(ncells_per_comm*sizeof(PetscScalar *));
  PetscScalar **vzcy_cell=(PetscScalar **)malloc(ncells_per_comm*sizeof(PetscScalar *));
  for(int i=0;i<ncells_per_comm;i++){
    bx_cell[i]=(PetscScalar *)malloc(ng2_per_cell*sizeof(PetscScalar));
    by_cell[i]=(PetscScalar *)malloc(ng2_per_cell*sizeof(PetscScalar));
    vxcx_cell[i]=(PetscScalar *)malloc(ng2_per_cell*sizeof(PetscScalar));
    vycx_cell[i]=(PetscScalar *)malloc(ng2_per_cell*sizeof(PetscScalar));
    vzcx_cell[i]=(PetscScalar *)malloc(ng2_per_cell*sizeof(PetscScalar));
    vxcy_cell[i]=(PetscScalar *)malloc(ng2_per_cell*sizeof(PetscScalar));
    vycy_cell[i]=(PetscScalar *)malloc(ng2_per_cell*sizeof(PetscScalar));
    vzcy_cell[i]=(PetscScalar *)malloc(ng2_per_cell*sizeof(PetscScalar));
  }

  double Lx_global=Lx*ncells_along_x, Ly_global=Ly*ncells_along_y;
  create_incident_pw(bx_natural, by_natural, freq, n_medium, polar,azimuth, input_polarization, Lx_global,Ly_global, Ngrid);
  MPI_Barrier(PETSC_COMM_WORLD);
  
  double Lxy[2]={Lx_global,Ly_global};
  MPI_Barrier(PETSC_COMM_WORLD);
  create_near2far(vxcx_natural, vycx_natural, vzcx_natural,
		  vxcy_natural, vycy_natural, vzcy_natural,
		  Ngrid, Lxy, oxy, symxy,
		  xyzfar,
		  freq, 1.0,1.0);
  MPI_Barrier(PETSC_COMM_WORLD);
  
  tandem2natural(bx_tandem, bx_natural, ngrid[0], ngrid[1], 1, ncells_along_x, ncells_along_y, -1);
  tandem2comm(colour, ng2_per_comm, bx_tandem, bx_comm);
  comm2cells(ng2_per_cell, ncells_per_comm, bx_comm, bx_cell,  1);
  tandem2natural(by_tandem, by_natural, ngrid[0], ngrid[1], 1, ncells_along_x, ncells_along_y, -1);
  tandem2comm(colour, ng2_per_comm, by_tandem, by_comm);
  comm2cells(ng2_per_cell, ncells_per_comm, by_comm, by_cell,  1);

  tandem2natural(vxcx_tandem, vxcx_natural, ngrid[0], ngrid[1], 1, ncells_along_x, ncells_along_y, -1);
  tandem2comm(colour, ng2_per_comm, vxcx_tandem, vxcx_comm);
  comm2cells(ng2_per_cell, ncells_per_comm, vxcx_comm, vxcx_cell,  1);
  tandem2natural(vycx_tandem, vycx_natural, ngrid[0], ngrid[1], 1, ncells_along_x, ncells_along_y, -1);
  tandem2comm(colour, ng2_per_comm, vycx_tandem, vycx_comm);
  comm2cells(ng2_per_cell, ncells_per_comm, vycx_comm, vycx_cell,  1);
  tandem2natural(vzcx_tandem, vzcx_natural, ngrid[0], ngrid[1], 1, ncells_along_x, ncells_along_y, -1);
  tandem2comm(colour, ng2_per_comm, vzcx_tandem, vzcx_comm);
  comm2cells(ng2_per_cell, ncells_per_comm, vzcx_comm, vzcx_cell,  1);

  tandem2natural(vxcy_tandem, vxcy_natural, ngrid[0], ngrid[1], 1, ncells_along_x, ncells_along_y, -1);
  tandem2comm(colour, ng2_per_comm, vxcy_tandem, vxcy_comm);
  comm2cells(ng2_per_cell, ncells_per_comm, vxcy_comm, vxcy_cell,  1);
  tandem2natural(vycy_tandem, vycy_natural, ngrid[0], ngrid[1], 1, ncells_along_x, ncells_along_y, -1);
  tandem2comm(colour, ng2_per_comm, vycy_tandem, vycy_comm);
  comm2cells(ng2_per_cell, ncells_per_comm, vycy_comm, vycy_cell,  1);
  tandem2natural(vzcy_tandem, vzcy_natural, ngrid[0], ngrid[1], 1, ncells_along_x, ncells_along_y, -1);
  tandem2comm(colour, ng2_per_comm, vzcy_tandem, vzcy_comm);
  comm2cells(ng2_per_cell, ncells_per_comm, vzcy_comm, vzcy_cell,  1);

  for(int i=0;i<ncells_per_comm;i++){
    for(int j=0;j<ng2_per_cell;j++){
      cell[i].bx[j]=bx_cell[i][j];
      cell[i].by[j]=by_cell[i][j];
      cell[i].vx_cx[j]=vxcx_cell[i][j];
      cell[i].vy_cx[j]=vycx_cell[i][j];
      cell[i].vz_cx[j]=vzcx_cell[i][j];
      cell[i].vx_cy[j]=vxcy_cell[i][j];
      cell[i].vy_cy[j]=vycy_cell[i][j];
      cell[i].vz_cy[j]=vzcy_cell[i][j];
    }
  }

  
  for(int i=0;i<ncells_per_comm;i++){
    free(bx_cell[i]);
    free(by_cell[i]);
    free(vxcx_cell[i]);
    free(vycx_cell[i]);
    free(vzcx_cell[i]);
    free(vxcy_cell[i]);
    free(vycy_cell[i]);
    free(vzcy_cell[i]);
  }
  free(bx_cell);
  free(by_cell);
  free(vxcx_cell);
  free(vycx_cell);
  free(vzcx_cell);
  free(vxcy_cell);
  free(vycy_cell);
  free(vzcy_cell);
    
  free(bx_natural);
  free(bx_tandem);
  free(bx_comm);
  free(by_natural);
  free(by_tandem);
  free(by_comm);
  free(vxcx_natural);
  free(vxcx_tandem);
  free(vxcx_comm);
  free(vycx_natural);
  free(vycx_tandem);
  free(vycx_comm);
  free(vzcx_natural);
  free(vzcx_tandem);
  free(vzcx_comm);
  free(vxcy_natural);
  free(vxcy_tandem);
  free(vxcy_comm);
  free(vycy_natural);
  free(vycy_tandem);
  free(vycy_comm);
  free(vzcy_natural);
  free(vzcy_tandem);
  free(vzcy_comm);

}


