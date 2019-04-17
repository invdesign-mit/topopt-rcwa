#include <complex>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "near2far.h"

using namespace std;

void green3d(std::complex<double> *EH, const double *x,
	     double freq, double eps, double mu,
	     const double *x0, int comp, std::complex<double> f0)
{

  double rhat[3] = {x[0]-x0[0],x[1]-x0[1],x[2]-x0[2]};
  double r = sqrt(rhat[0]*rhat[0]+rhat[1]*rhat[1]+rhat[2]*rhat[2]);
  rhat[0]=rhat[0]/r, rhat[1]=rhat[1]/r, rhat[2]=rhat[2]/r;

  double n = sqrt(eps*mu);
  double k = 2*M_PI*freq*n;
  std::complex<double> ikr = std::complex<double>(0.0, k*r);
  double ikr2   = -(k*r)*(k*r);
  /* note that SCUFF-EM computes the fields from the dipole moment p,
       whereas we need it from the current J = -i*omega*p, so our result
       is divided by -i*omega compared to SCUFF */
  std::complex<double> expfac = f0 * (k*n/(4*M_PI*r)) * exp(std::complex<double>(0,k*r + M_PI*0.5));
  double Z = sqrt(mu/eps);

  double p[3]={0,0,0};
  p[comp%3]=1.0;
  double pdotrhat = p[0]*rhat[0] + p[1]*rhat[1] + p[2]*rhat[2];
  
  double rhatcrossp[3] = {rhat[1] * p[2] -
			  rhat[2] * p[1],
			  rhat[2] * p[0] -
			  rhat[0] * p[2],
			  rhat[0] * p[1] -
			  rhat[1] * p[0]};
  
  /* compute the various scalar quantities in the point source formulae */
  std::complex<double> term1 =  1.0 - 1.0/ikr + 1.0/ikr2;
  std::complex<double> term2 = (-1.0 + 3.0/ikr - 3.0/ikr2) * pdotrhat;
  std::complex<double> term3 = (1.0 - 1.0/ikr);
  /* now assemble everything based on source type */
  if (comp<3) {
    expfac /= eps;
    EH[0] = expfac * (term1*p[0] + term2*rhat[0]);
    EH[1] = expfac * (term1*p[1] + term2*rhat[1]);
    EH[2] = expfac * (term1*p[2] + term2*rhat[2]);
    EH[3] = expfac*term3*rhatcrossp[0] / Z;
    EH[4] = expfac*term3*rhatcrossp[1] / Z;
    EH[5] = expfac*term3*rhatcrossp[2] / Z;
  }
  else {
    expfac /= mu;
    EH[0] = -expfac*term3*rhatcrossp[0] * Z;
    EH[1] = -expfac*term3*rhatcrossp[1] * Z;
    EH[2] = -expfac*term3*rhatcrossp[2] * Z;
    EH[3] = expfac * (term1*p[0] + term2*rhat[0]);
    EH[4] = expfac * (term1*p[1] + term2*rhat[1]);
    EH[5] = expfac * (term1*p[2] + term2*rhat[2]);
  }

  //conjugation
  EH[0]=conj(EH[0]);
  EH[1]=conj(EH[1]);
  EH[2]=conj(EH[2]);
  EH[3]=conj(EH[3]);
  EH[4]=conj(EH[4]);
  EH[5]=conj(EH[5]);
  
}


void create_near2far(PetscScalar *vx_cx, PetscScalar *vy_cx, PetscScalar *vz_cx,
		     PetscScalar *vx_cy, PetscScalar *vy_cy, PetscScalar *vz_cy,
		     int nxy[2], double Lxy[2], double oxy[2], int symxy[2],
		     double xyzfar[3],
		     double freq, double eps, double mu)
{

  double dg[2]={Lxy[0]/nxy[0],Lxy[1]/nxy[1]};
  PetscScalar Fex[6],Fey[6];


  for(int iy=0;iy<nxy[1];iy++){
    for(int ix=0;ix<nxy[0];ix++){

      int i=ix+nxy[0]*iy;
      PetscScalar Mx =  1.0; //get the right magnetic currents
      PetscScalar My = -1.0; //get the right magnetic currents

      double xyznear[3]={ (ix+0.5)*dg[0] - oxy[0], (iy+0.5)*dg[1] - oxy[1], 0.0 };
      green3d(Fey, xyzfar, freq,eps,mu, xyznear, 3, Mx);
      green3d(Fex, xyzfar, freq,eps,mu, xyznear, 4, My);
      vx_cx[i]=Fex[0], vx_cy[i]=Fey[0];
      vy_cx[i]=Fex[1], vy_cy[i]=Fey[1];
      vz_cx[i]=Fex[2], vz_cy[i]=Fey[2];

      if(symxy[0]==1){ //mirror fields at (-x,y,z); even boundary 
	double xyz[3]={ -xyznear[0], xyznear[1], xyznear[2] };
	green3d(Fey, xyzfar, freq,eps,mu, xyz, 3, Mx);
	green3d(Fex, xyzfar, freq,eps,mu, xyz, 4,-My);
	vx_cx[i]+=Fex[0], vx_cy[i]+=Fey[0];
	vy_cx[i]+=Fex[1], vy_cy[i]+=Fey[1];
	vz_cx[i]+=Fex[2], vz_cy[i]+=Fey[2];
      }
      if(symxy[0]==-1){ //mirror fields at (-x,y,z); odd boundary 
	double xyz[3]={ -xyznear[0], xyznear[1], xyznear[2] };
	green3d(Fey, xyzfar, freq,eps,mu, xyz, 3,-Mx);
	green3d(Fex, xyzfar, freq,eps,mu, xyz, 4, My);
	vx_cx[i]+=Fex[0], vx_cy[i]+=Fey[0];
	vy_cx[i]+=Fex[1], vy_cy[i]+=Fey[1];
	vz_cx[i]+=Fex[2], vz_cy[i]+=Fey[2];
      }
      if(symxy[1]==1){ //mirror fields at (x,-y,z); even boundary 
	double xyz[3]={  xyznear[0],-xyznear[1], xyznear[2] };
	green3d(Fey, xyzfar, freq,eps,mu, xyz, 3,-Mx);
	green3d(Fex, xyzfar, freq,eps,mu, xyz, 4, My);
	vx_cx[i]+=Fex[0], vx_cy[i]+=Fey[0];
	vy_cx[i]+=Fex[1], vy_cy[i]+=Fey[1];
	vz_cx[i]+=Fex[2], vz_cy[i]+=Fey[2];
      }
      if(symxy[1]==-1){ //mirror fields at (x,-y,z); odd boundary 
	double xyz[3]={  xyznear[0],-xyznear[1], xyznear[2] };
	green3d(Fey, xyzfar, freq,eps,mu, xyz, 3, Mx);
	green3d(Fex, xyzfar, freq,eps,mu, xyz, 4,-My);
	vx_cx[i]+=Fex[0], vx_cy[i]+=Fey[0];
	vy_cx[i]+=Fex[1], vy_cy[i]+=Fey[1];
	vz_cx[i]+=Fex[2], vz_cy[i]+=Fey[2];
      }
      if( (symxy[0]==1 && symxy[1]==1) || (symxy[0]==-1 && symxy[1]==-1) ){ //mirror fields at (-x,-y,z); only even/odd boundaries 
	double xyz[3]={ -xyznear[0],-xyznear[1], xyznear[2] };
	green3d(Fey, xyzfar, freq,eps,mu, xyz, 3,-Mx);
	green3d(Fex, xyzfar, freq,eps,mu, xyz, 4,-My);
	vx_cx[i]+=Fex[0], vx_cy[i]+=Fey[0];
	vy_cx[i]+=Fex[1], vy_cy[i]+=Fey[1];
	vz_cx[i]+=Fex[2], vz_cy[i]+=Fey[2];
      }
      if( (symxy[0]==1 && symxy[1]==-1) || (symxy[0]==-1 && symxy[1]==1) ){ //mirror fields at (-x,-y,z); mixed boundaries 
	double xyz[3]={ -xyznear[0],-xyznear[1], xyznear[2] };
	green3d(Fey, xyzfar, freq,eps,mu, xyz, 3, Mx);
	green3d(Fex, xyzfar, freq,eps,mu, xyz, 4, My);
	vx_cx[i]+=Fex[0], vx_cy[i]+=Fey[0];
	vy_cx[i]+=Fex[1], vy_cy[i]+=Fey[1];
	vz_cx[i]+=Fex[2], vz_cy[i]+=Fey[2];
      }

    }
  }
}




//**********************
void create_near2far_2d(PetscScalar *vx_cx, PetscScalar *vy_cx, PetscScalar *vz_cx,
			PetscScalar *vx_cy, PetscScalar *vy_cy, PetscScalar *vz_cy,
			int nxy[2], double Lxy[2], double oxy[2], int symxy[2],
			double xyzfar_3d[3],
			double freq, double eps, double mu)
{

  double dg[2]={Lxy[0]/nxy[0],Lxy[1]/nxy[1]};
  PetscScalar Fex[6];

  double xyzfar[3]={xyzfar_3d[0],xyzfar_3d[1],xyzfar_3d[2]};

  for(int iy=0;iy<nxy[1];iy++){
    for(int ix=0;ix<nxy[0];ix++){

      int i=ix+nxy[0]*iy;
      PetscScalar My = -1.0; //get the right magnetic currents

      double xyznear[3]={ (ix+0.5)*dg[0] - oxy[0], (iy+0.5)*dg[1] - oxy[1], 0.0 };
      xyzfar[0]=xyznear[0];
      green3d(Fex, xyzfar, freq,eps,mu, xyznear, 4, My);
      vx_cx[i]=Fex[0], vx_cy[i]=0.0;
      vy_cx[i]=0.0,    vy_cy[i]=0.0;
      vz_cx[i]=0.0,    vz_cy[i]=0.0;
    }
  }

}


void create_near2far_2d_phase(PetscScalar *vx_cx, PetscScalar *vy_cx, PetscScalar *vz_cx,
			      PetscScalar *vx_cy, PetscScalar *vy_cy, PetscScalar *vz_cy,
			      int nxy[2], double Lxy[2], double oxy[2], int symxy[2],
			      double xyzfar[3],
			      double freq, double eps, double mu)
{

  double dg[2]={Lxy[0]/nxy[0],Lxy[1]/nxy[1]};

  double yF=xyzfar[1];
  double zF=xyzfar[2];
  double n=sqrt(eps*mu);
  
  for(int iy=0;iy<nxy[1];iy++){
    for(int ix=0;ix<nxy[0];ix++){

      int i=ix+nxy[0]*iy;

      double xN=(ix+0.5)*dg[0]-oxy[0];
      double yN=(iy+0.5)*dg[1]-oxy[1];
      double zN=0.0;
      double xF=xN;

      double r=sqrt( pow(xF-xN,2) + pow(yF-yN,2) + pow(zF-zN,2) );
      double phi=2*M_PI*n*freq*(r-zF);
      vx_cx[i]=exp( -1.0 * std::complex<double>(0.0,phi) );

      vx_cy[i]=0.0;
      vy_cx[i]=0.0,    vy_cy[i]=0.0;
      vz_cx[i]=0.0,    vz_cy[i]=0.0;
    }
  }

}
