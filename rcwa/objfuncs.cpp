# include "rcwa.h"
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex>
# include "gsel.h"
# include "fft_iface.h"
# include <iostream>
# include <iomanip>
# include <cmath>
# include <string>
# include "numalloc.h"
# include <TBLAS.h>
#ifdef HAVE_BLAS
# include <TBLAS_ext.h>
#endif
# include "wrap.h"
# include "objfuncs.h"

using namespace std;
  
void private_grad(double *epsbar, double kbloch[2], rcwa_ *rcwa, cell_ *cell,
		  std::complex<double> *cx, std::complex<double> *cy,
		  std::complex<double> *grad_complex, std::complex<double> *intface_grad){

  int ng2=rcwa->ngrid[0]*rcwa->ngrid[1];
  double dg2=rcwa->dg[0]*rcwa->dg[1];
  
  rhs2a(rcwa->nG,rcwa->G, rcwa->ngrid, cx,cy, -kbloch[0],-kbloch[1], rcwa->omega, rcwa->dg, -1, rcwa->bNx_adj,rcwa->bNy_adj, 1);
  for(int i=0;i<rcwa->nlayers_patterned;i++){
    rcwa_solve_layer(rcwa,rcwa->layerID_patterned[i]);
    double t=rcwa->thickness[rcwa->layerID_patterned[i]];
    int nz=cell->nz_integrate_patterned[i];
    double dz=t/double(nz);
    if(intface_grad){ intface_grad[2*i+0]=0.0, intface_grad[2*i+1]=0.0; };
    for(int iz=0;iz<=nz;iz++){
      if(iz==0)
	rcwa_translate_and_getfields(rcwa,rcwa->layerID_patterned[i], 0,  t, kbloch);
      else
	rcwa_translate_and_getfields(rcwa,rcwa->layerID_patterned[i],dz,-dz, kbloch);
      for(int j=0;j<ng2;j++){
	std::complex<double> tmp=pow(rcwa->omega,2)*cell->epsdiff_patterned[i]*( rcwa->efield_adj[3*j+0]*rcwa->efield_fwd[3*j+0] + \
										 rcwa->efield_adj[3*j+1]*rcwa->efield_fwd[3*j+1] + \
										 rcwa->efield_adj[3*j+2]*rcwa->efield_fwd[3*j+2] );
	if(iz==0)
	  grad_complex[j+i*ng2]  = dz/3.0 * tmp * dg2;
	else if(iz>0 && iz<nz)
	  grad_complex[j+i*ng2] += dz/3.0 * ((iz%2)==0 ? 2.0 : 4.0) * tmp * dg2;
	else
	  grad_complex[j+i*ng2] += dz/3.0 * tmp * dg2;
	
	if(intface_grad){
	  if(iz==0){
	    double eps1=cell->epsdiff_patterned[i]*epsbar[j+i*ng2]+cell->epsbkg_patterned[i];
	    double eps2=cell->n_uniform[i]*cell->n_uniform[i];
	    intface_grad[2*i+0]+=std::pow(rcwa->omega,2)*( (eps2-eps1) * rcwa->efield_adj[3*j+0] * rcwa->efield_fwd[3*j+0] + \
							   (eps2-eps1) * rcwa->efield_adj[3*j+1] * rcwa->efield_fwd[3*j+1] ) * dg2; 
	    //(1.0/eps1 - 1.0/eps2) * Q.dfield[3*j+2] * F.dfield[3*j+2] ) * dg2;
	  }
	  if(iz==nz){
	    double eps1=cell->n_uniform[i+1]*cell->n_uniform[i+1];
	    double eps2=cell->epsdiff_patterned[i]*epsbar[j+i*ng2]+cell->epsbkg_patterned[i];
	    intface_grad[2*i+1]+=std::pow(rcwa->omega,2)*( (eps2-eps1) * rcwa->efield_adj[3*j+0] * rcwa->efield_fwd[3*j+0] + \
							   (eps2-eps1) * rcwa->efield_adj[3*j+1] * rcwa->efield_fwd[3*j+1] ) * dg2; 
	    //(1.0/eps1 - 1.0/eps2) * Q.dfield[3*j+2] * F.dfield[3*j+2] ) * dg2;
	  }
	}

      }
    }
  }

}

void v_dot_E(double *epsbar, cell_ *cell,
	     std::complex<double> *bx, std::complex<double> *by,
	     std::complex<double> *ux, std::complex<double> *uy, double uobj[2], std::complex<double> *ugrad, std::complex<double> *intface_ugrad,
	     std::complex<double> *vx, std::complex<double> *vy, double vobj[2], std::complex<double> *vgrad, std::complex<double> *intface_vgrad,
	     std::complex<double> *wx, std::complex<double> *wy, double wobj[2], std::complex<double> *wgrad, std::complex<double> *intface_wgrad){

  rcwa_ rcwa;
  rcwa.nlayers_uniform=cell->nlayers_uniform;
  rcwa.nlayers_patterned=cell->nlayers_patterned;
  rcwa.nG=cell->nG;
  rcwa.Lr[0]=cell->Lr[0], rcwa.Lr[1]=cell->Lr[1], rcwa.Lr[2]=cell->Lr[2], rcwa.Lr[3]=cell->Lr[3];
  rcwa.ngrid[0]=cell->ngrid[0];
  rcwa.ngrid[1]=cell->ngrid[1];
  rcwa.freq=cell->freq;
  double kbloch[2];
  kbloch[0]=2 * M_PI * cell->freq * cell->dir_cosine[0];
  kbloch[1]=2 * M_PI * cell->freq * cell->dir_cosine[1];
  
  memset_rcwa(&rcwa,epsbar,kbloch,
	      cell->layerID_uniform,cell->n_uniform,
	      cell->layerID_patterned,cell->epsdiff_patterned,cell->epsbkg_patterned,
	      cell->thickness);

  int ng2=rcwa.ngrid[0]*rcwa.ngrid[1];
  double dg2=rcwa.dg[0]*rcwa.dg[1];
  rhs2a(rcwa.nG,rcwa.G, rcwa.ngrid, bx,by, kbloch[0],kbloch[1], rcwa.omega, rcwa.dg, 1, rcwa.a0x_fwd,rcwa.a0y_fwd, 1);
  rcwa_getfields_fwd(&rcwa, rcwa.nlayers-1, rcwa.thickness[rcwa.nlayers-1], kbloch);
  std::complex<double>tmpobju(0,0);
  std::complex<double>tmpobjv(0,0);
  std::complex<double>tmpobjw(0,0);
  for(int i=0;i<ng2;i++){
    if(ux && uy) tmpobju += dg2 * (ux[i]*rcwa.efield_fwd[3*i+0] + uy[i]*rcwa.efield_fwd[3*i+1]);
    if(vx && vy) tmpobjv += dg2 * (vx[i]*rcwa.efield_fwd[3*i+0] + vy[i]*rcwa.efield_fwd[3*i+1]);
    if(wx && wy) tmpobjw += dg2 * (wx[i]*rcwa.efield_fwd[3*i+0] + wy[i]*rcwa.efield_fwd[3*i+1]);
  }
  uobj[0]=real(tmpobju), uobj[1]=imag(tmpobju);
  vobj[0]=real(tmpobjv), vobj[1]=imag(tmpobjv);
  wobj[0]=real(tmpobjw), wobj[1]=imag(tmpobjw);

  if(ux && uy) private_grad(epsbar,kbloch,&rcwa,cell,ux,uy,ugrad,intface_ugrad);
  if(vx && vy) private_grad(epsbar,kbloch,&rcwa,cell,vx,vy,vgrad,intface_vgrad);
  if(wx && wy) private_grad(epsbar,kbloch,&rcwa,cell,wx,wy,wgrad,intface_wgrad);

  destroy_rcwa(&rcwa);
  
}



