/* Copyright (C) 2009-2011, Stanford University
 * This file is part of S4
 * Written by Victor Liu (vkl@stanford.edu)
 *
 * S4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * S4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "fft_iface.h"
#include <cstdlib>

#include <fftw3.h>

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
static pthread_mutex_t mutex;
#endif

#include <cmath>
#include <TBLAS.h>
#ifdef HAVE_BLAS
# include <TBLAS_ext.h>
#endif
#include <LinearSolve.h>
#ifdef HAVE_LAPACK
# include <LinearSolve_lapack.h>
# include <Eigensystems_lapack.h>
#else
# include <Eigensystems.h>
#endif
#include <limits>
#include "numalloc.h"

# include <iostream>
# include <iomanip>

int fft_next_fast_size(int n){
    while(1){
        int m=n;
        while( (m%2) == 0 ) m/=2;
        while( (m%3) == 0 ) m/=3;
        while( (m%5) == 0 ) m/=5;
        if(m<=1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}

std::complex<double> *fft_alloc_complex(size_t n){
  return (std::complex<double>*)(fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n);
}

void fft_free(void *p){
  fftw_free(p);
}

struct tag_fft_plan{
  fftw_plan plan;
};

fft_plan fft_plan_dft_2d(
	int n[2],
	std::complex<double> *in, std::complex<double> *out,
	int sign
){
  fft_plan plan = NULL;
# ifdef HAVE_LIBPTHREAD
  pthread_mutex_lock(&mutex);
# endif
  fftw_plan p;
  p = fftw_plan_dft(2, n, (fftw_complex*)in, (fftw_complex*)out, sign, FFTW_ESTIMATE);
# ifdef HAVE_LIBPTHREAD
  pthread_mutex_unlock(&mutex);
# endif
  if(NULL != p){
    plan = (fft_plan)malloc(sizeof(tag_fft_plan));
    plan->plan = p;
  }
  return plan;
}

void fft_plan_exec(const fft_plan plan){
  fftw_execute(plan->plan);
}

void fft_plan_destroy(fft_plan plan){
  if(NULL == plan){ return; }
# ifdef HAVE_LIBPTHREAD
  pthread_mutex_lock(&mutex);
# endif
  fftw_destroy_plan(plan->plan);
# ifdef HAVE_LIBPTHREAD
  pthread_mutex_unlock(&mutex);
# endif
  free(plan);
}

void fft_init(){
#ifdef HAVE_LIBPTHREAD
  if(pthread_mutex_init(&mutex, NULL)){
    printf("Unable to initialize a mutex for FFT module\n");
  }
#endif
}

void fft_destroy(){
  fftw_cleanup();
#ifdef HAVE_LIBPTHREAD
  pthread_mutex_destroy(&mutex);
#endif
}

int FMMGetEpsilon_FFT(const int n, int *ngrid, const std::complex<double> *eps, const int *G, std::complex<double> *Epsilon2, std::complex<double> *Epsilon_inv){

  const int n2 = 2*n;
  
  const int ng2 = ngrid[0]*ngrid[1];
  const double ing2 = 1./ng2;
  // The grid needs to hold 5 matrix elements: xx,xy,yx,yy,zz
  // We actually make 5 different grids to facilitate the fft routines

  std::complex<double> *work = (std::complex<double>*)my_malloc(sizeof(std::complex<double>)*(6*ng2));
  std::complex<double>*fxx = work;
  std::complex<double>*fxy = fxx + ng2;
  std::complex<double>*fyx = fxy + ng2;
  std::complex<double>*fyy = fyx + ng2;
  std::complex<double>*fzz = fyy + ng2;
  std::complex<double>*Fto = fzz + ng2;

  for(int i = 0; i < ng2; i++){
    fxx[i]=eps[i];
    fyy[i]=eps[i];
    fzz[i]=eps[i];
    fxy[i]=0;
    fyx[i]=0;
  }

  int ngrid_fftw[2]={ngrid[1],ngrid[0]}; //nx = ngrid[0] is the faster growing index which must be placed later to match with fftw.
  fft_plan plans[5];
  for(int i = 0; i <= 4; ++i){
    plans[i] = fft_plan_dft_2d(ngrid_fftw, fxx+i*ng2, Fto, 1);
  }

  fft_plan_exec(plans[4]);
  for(int j = 0; j < n; ++j){
    for(int i = 0; i < n; ++i){
      int f[2] = {G[2*i+0]-G[2*j+0],G[2*i+1]-G[2*j+1]};
      if(f[0] < 0){ f[0] += ngrid[0]; }
      if(f[1] < 0){ f[1] += ngrid[1]; }
      Epsilon2[i+j*n] = ing2 * Fto[f[0]+f[1]*ngrid[0]];
    }
  }

  // Epsilon_inv needs inverting
  RNP::TBLAS::SetMatrix<'A'>(n,n, 0.,1., Epsilon_inv,n);
  int solve_info;
  RNP::LinearSolve<'N'>(n,n, Epsilon2,n, Epsilon_inv,n, &solve_info, NULL);

  // We fill in the quarter blocks of F in Fortran order
  for(int w = 0; w < 4; ++w){
    int Ecol = (w&1 ? n : 0);
    int Erow = (w&2 ? n : 0);


    fft_plan_exec(plans[w]);

    for(int j = 0; j < n; ++j){
      for(int i = 0; i < n; ++i){
	int f[2] = {G[2*i+0]-G[2*j+0],G[2*i+1]-G[2*j+1]};
	if(f[0] < 0){ f[0] += ngrid[0]; }
	if(f[1] < 0){ f[1] += ngrid[1]; }
	Epsilon2[Erow+i+(Ecol+j)*n2] = ing2 * Fto[f[0]+f[1]*ngrid[0]];
      }
    }
  }
  for(int i = 0; i <= 4; ++i){
    fft_plan_destroy(plans[i]);
  }

  my_free(work);
  return 0;
}

//FDFD assumes e^[ i w t] ==> so forward wave is e^[-i k.r + i w t]
//RCWA assumes e^[-i w t] ==> so forward wave is e^[ i k.r - i w t]
//RCWA assumes e^[-i w t] ==> curl curl E - omega^2 eps E = i omega J = rhs ==> J = rhs/(i omega) = conj[rhs_fdfd]/(i omega)
int rhs2a(const int n, int *G, int *ngrid, std::complex<double> *rhsfdfd_x, std:: complex<double> *rhsfdfd_y, const double kx, const double ky, const double omega_real, const double *dg, const std::complex<double> extra_factor, std::complex<double> *ax, std::complex<double> *ay, int plan_dir){

  const int ng2=ngrid[0]*ngrid[1];
  const double ing2=1.0/ng2;

  const std::complex<double> imagI(0.0,1.0);

  //remove the Bloch factor exp(i * k_bloch \dot (x,y) ) from J
  std::complex<double> *Jx = (std::complex<double>*)malloc(sizeof(std::complex<double>)*ng2);
  std::complex<double> *Jy = (std::complex<double>*)malloc(sizeof(std::complex<double>)*ng2);
  for(int j=0;j<ngrid[1];j++){
    for(int i=0;i<ngrid[0];i++){
      int ind = i+j*ngrid[0];
      double x = (i+0.5)*dg[0];
      double y = (j+0.5)*dg[1];
      std::complex<double> rj=std::exp( imagI * (kx*x + ky*y) ) * imagI*omega_real;
      Jx[ind] = std::conj(rhsfdfd_x[ind])/rj;
      Jy[ind] = std::conj(rhsfdfd_y[ind])/rj;

    }
  }

  int ngrid_fftw[2]={ngrid[1],ngrid[0]}; //nx=ngrid[0] is the faster growing index which must be placed later to match with fftw.
  std::complex<double> *jx = (std::complex<double>*)malloc(sizeof(std::complex<double>)*ng2);
  std::complex<double> *jy = (std::complex<double>*)malloc(sizeof(std::complex<double>)*ng2);
  fft_plan plan_x, plan_y;
  plan_x = fft_plan_dft_2d(ngrid_fftw, Jx, jx, plan_dir); //sign_debug plan_dir=1 for b, -1 for v
  plan_y = fft_plan_dft_2d(ngrid_fftw, Jy, jy, plan_dir); //sign_debug plan_dir=1 for b, -1 for v
  fft_plan_exec(plan_x);
  fft_plan_exec(plan_y);

  for(int i=0;i<n;i++){
    int f[2]={G[2*i+0],G[2*i+1]};
    if(f[0] < 0){ f[0] += ngrid[0]; }
    if(f[1] < 0){ f[1] += ngrid[1]; }
    ax[i] = ( jy[f[0] + f[1]*ngrid[0]]/2.0) * extra_factor * ing2;
    ay[i] = (-jx[f[0] + f[1]*ngrid[0]]/2.0) * extra_factor * ing2;
  }

  fft_plan_destroy(plan_x);
  fft_plan_destroy(plan_y);
  free(Jx);
  free(Jy);
  free(jx);
  free(jy);
  return 0;

}

