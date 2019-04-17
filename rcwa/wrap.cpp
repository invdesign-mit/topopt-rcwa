# include "rcwa.h"
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <complex>
# include "gsel.h"
# include "fft_iface.h"
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <string>
# include "numalloc.h"
# include <TBLAS.h>
#ifdef HAVE_BLAS
# include <TBLAS_ext.h>
#endif
# include "wrap.h"

void memset_rcwa(rcwa_ *rcwa, double *epsbar, double kbloch[2],
		 int *layerID_uniform, double *n_uniform,
		 int *layerID_patterned, double *epsdiff, double *epsbkg,
		 double *thickness){

  rcwa->kbloch[0] = -kbloch[0];
  rcwa->kbloch[1] = -kbloch[1];
  
  Lattice_Reciprocate(rcwa->Lr, rcwa->Lk);
  rcwa->G = (int *)my_malloc(2*rcwa->nG*sizeof(int));
  G_select(0,&(rcwa->nG),rcwa->Lk,rcwa->G);
  rcwa->omega=2.0*M_PI*rcwa->freq;
  rcwa->omega_complex=std::complex<double>(rcwa->omega,0.0);
  rcwa->kx_fwd=(double *)malloc(rcwa->nG*sizeof(double));
  rcwa->ky_fwd=(double *)malloc(rcwa->nG*sizeof(double));
  rcwa->kx_adj=(double *)malloc(rcwa->nG*sizeof(double));
  rcwa->ky_adj=(double *)malloc(rcwa->nG*sizeof(double));
  Lattice_SetKs(rcwa->nG, rcwa->kbloch[0], rcwa->kbloch[1], rcwa->G, rcwa->Lk, rcwa->kx_fwd, rcwa->ky_fwd);
  Lattice_SetKs(rcwa->nG,-rcwa->kbloch[0],-rcwa->kbloch[1], rcwa->G, rcwa->Lk, rcwa->kx_adj, rcwa->ky_adj);
  
  rcwa->nlayers=rcwa->nlayers_uniform+rcwa->nlayers_patterned;
  rcwa->layerID_uniform=(int *) malloc(rcwa->nlayers_uniform * sizeof(int));
  rcwa->layerID_patterned=(int *) malloc(rcwa->nlayers_patterned * sizeof(int));
  rcwa->thickness=(double *)malloc(rcwa->nlayers * sizeof(double));
  for(int i=0;i<rcwa->nlayers_uniform;i++){
    rcwa->layerID_uniform[i]=layerID_uniform[i];
    rcwa->thickness[layerID_uniform[i]]=thickness[layerID_uniform[i]];
  }
  for(int i=0;i<rcwa->nlayers_patterned;i++){
    rcwa->layerID_patterned[i]=layerID_patterned[i];
    rcwa->thickness[layerID_patterned[i]]=thickness[layerID_patterned[i]];
  }
    
  int ng2=rcwa->ngrid[0]*rcwa->ngrid[1];
  rcwa->dg[0]=rcwa->Lr[0]/double(rcwa->ngrid[0]);
  rcwa->dg[1]=rcwa->Lr[3]/double(rcwa->ngrid[1]);
  
  rcwa->eps_uniform = (std::complex<double> *) malloc(ng2*sizeof(std::complex<double>));
  rcwa->eps_uniform_inv = (std::complex<double> *) malloc(ng2*sizeof(std::complex<double>));
  if(n_uniform!=NULL){
    for(int i=0;i<rcwa->nlayers_uniform;i++){
      rcwa->eps_uniform[i] = std::complex<double>(n_uniform[i]*n_uniform[i],0.0);
      rcwa->eps_uniform_inv[i] = 1.0/rcwa->eps_uniform[i];
    }
  }else{
    std::cout << "ERROR: n_uniform must not be NULL.\n";
  }

  rcwa->eps_patterned = (std::complex<double> **) malloc(rcwa->nlayers_patterned*sizeof(std::complex<double>*));
  for(int j=0;j<rcwa->nlayers_patterned;j++){
    rcwa->eps_patterned[j] = (std::complex<double> *) malloc(ng2*sizeof(std::complex<double>));
    if(epsbar!=NULL && epsdiff!=NULL && epsbkg!=NULL)
      for(int i=0;i<ng2;i++)
	rcwa->eps_patterned[j][i] = std::complex<double>(epsdiff[j]*epsbar[i+j*ng2]+epsbkg[j],0.0);
  }

  rcwa->Epsilon2=(std::complex<double>**)malloc(rcwa->nlayers_patterned*sizeof(std::complex<double>*));
  rcwa->Epsilon_inv=(std::complex<double>**)malloc(rcwa->nlayers_patterned*sizeof(std::complex<double>*));
  rcwa->epstype=(int *)malloc(rcwa->nlayers*sizeof(int));
  rcwa->eps_inv_list=(const std::complex<double>**)malloc(rcwa->nlayers*sizeof(std::complex<double>*));

  for(int i=0;i<rcwa->nlayers_patterned;i++){
    rcwa->Epsilon2[i] = (std::complex<double>*) my_malloc(2*rcwa->nG*2*rcwa->nG*sizeof(std::complex<double>));
    rcwa->Epsilon_inv[i] = (std::complex<double>*) my_malloc(rcwa->nG*rcwa->nG*sizeof(std::complex<double>));
    rcwa->eps_inv_list[rcwa->layerID_patterned[i]]=rcwa->Epsilon_inv[i];
    rcwa->epstype[rcwa->layerID_patterned[i]]=EPSILON2_TYPE_FULL;
  }
  for(int i=0;i<rcwa->nlayers_uniform;i++){
    rcwa->eps_inv_list[rcwa->layerID_uniform[i]]=&(rcwa->eps_uniform_inv[i]);
    rcwa->epstype[rcwa->layerID_uniform[i]]=EPSILON2_TYPE_BLKDIAG1_SCALAR;
  }

  rcwa->q_fwd   = (std::complex<double>**)malloc(rcwa->nlayers*sizeof(std::complex<double>*));
  rcwa->kp_fwd  = (std::complex<double>**)malloc(rcwa->nlayers*sizeof(std::complex<double>*));
  rcwa->phi_fwd = (std::complex<double>**)malloc(rcwa->nlayers*sizeof(std::complex<double>*));
  rcwa->q_adj   = (std::complex<double>**)malloc(rcwa->nlayers*sizeof(std::complex<double>*));
  rcwa->kp_adj  = (std::complex<double>**)malloc(rcwa->nlayers*sizeof(std::complex<double>*));
  rcwa->phi_adj = (std::complex<double>**)malloc(rcwa->nlayers*sizeof(std::complex<double>*));
  for(int i=0;i<rcwa->nlayers;i++){
    rcwa->q_fwd[i] = (std::complex<double>*)my_malloc(2*rcwa->nG*sizeof(std::complex<double>));
    rcwa->kp_fwd[i] = (std::complex<double>*)my_malloc(4*rcwa->nG*rcwa->nG*sizeof(std::complex<double>));
    rcwa->phi_fwd[i] = (std::complex<double>*)my_malloc(4*rcwa->nG*rcwa->nG*sizeof(std::complex<double>));
    rcwa->q_adj[i] = (std::complex<double>*)my_malloc(2*rcwa->nG*sizeof(std::complex<double>));
    rcwa->kp_adj[i] = (std::complex<double>*)my_malloc(4*rcwa->nG*rcwa->nG*sizeof(std::complex<double>));
    rcwa->phi_adj[i] = (std::complex<double>*)my_malloc(4*rcwa->nG*rcwa->nG*sizeof(std::complex<double>));
  }
  for(int i=0;i<rcwa->nlayers_uniform;i++){

    SolveLayerEigensystem_uniform(
				  rcwa->omega_complex,
				  rcwa->nG,
				  rcwa->kx_fwd,
				  rcwa->ky_fwd,
				  rcwa->eps_uniform[i],
				  rcwa->q_fwd[rcwa->layerID_uniform[i]],
				  rcwa->kp_fwd[rcwa->layerID_uniform[i]],
				  rcwa->phi_fwd[rcwa->layerID_uniform[i]]
				  );
    SolveLayerEigensystem_uniform(
				  rcwa->omega_complex,
				  rcwa->nG,
				  rcwa->kx_adj,
				  rcwa->ky_adj,
				  rcwa->eps_uniform[i],
				  rcwa->q_adj[rcwa->layerID_uniform[i]],
				  rcwa->kp_adj[rcwa->layerID_uniform[i]],
				  rcwa->phi_adj[rcwa->layerID_uniform[i]]
				  );

  }

  rcwa->ab_fwd  = (std::complex<double>*)my_malloc(4*rcwa->nG*sizeof(std::complex<double>));
  rcwa->a0_fwd  = &(rcwa->ab_fwd[0*2*rcwa->nG]);
  rcwa->bN_fwd  = &(rcwa->ab_fwd[1*2*rcwa->nG]);
  rcwa->a0x_fwd = &(rcwa->ab_fwd[0*1*rcwa->nG]);
  rcwa->a0y_fwd = &(rcwa->ab_fwd[1*1*rcwa->nG]);
  rcwa->bNx_fwd = &(rcwa->ab_fwd[2*1*rcwa->nG]);
  rcwa->bNy_fwd = &(rcwa->ab_fwd[3*1*rcwa->nG]);
    
  rcwa->ab_adj  = (std::complex<double>*)my_malloc(4*rcwa->nG*sizeof(std::complex<double>));
  rcwa->a0_adj  = &(rcwa->ab_adj[0*2*rcwa->nG]);
  rcwa->bN_adj  = &(rcwa->ab_adj[1*2*rcwa->nG]);
  rcwa->a0x_adj = &(rcwa->ab_adj[0*1*rcwa->nG]);
  rcwa->a0y_adj = &(rcwa->ab_adj[1*1*rcwa->nG]);
  rcwa->bNx_adj = &(rcwa->ab_adj[2*1*rcwa->nG]);
  rcwa->bNy_adj = &(rcwa->ab_adj[3*1*rcwa->nG]);

  for(int i=0;i<2*int(rcwa->nG);i++){
    rcwa->bN_fwd[i]=0.0;
    rcwa->a0_adj[i]=0.0;
  }

  rcwa->ab_internal_fwd  = (std::complex<double>*)my_malloc(4*rcwa->nG*sizeof(std::complex<double>));
  rcwa->ab_internal_adj  = (std::complex<double>*)my_malloc(4*rcwa->nG*sizeof(std::complex<double>));
  
  rcwa->efield_fwd=(std::complex<double>*)my_malloc(3*ng2*sizeof(std::complex<double>));
  rcwa->hfield_fwd=(std::complex<double>*)my_malloc(3*ng2*sizeof(std::complex<double>));
  rcwa->efield_adj=(std::complex<double>*)my_malloc(3*ng2*sizeof(std::complex<double>));
  rcwa->hfield_adj=(std::complex<double>*)my_malloc(3*ng2*sizeof(std::complex<double>));

  for(int i=0;i<rcwa->nlayers_patterned;i++){
    FMMGetEpsilon_FFT(rcwa->nG, rcwa->ngrid, rcwa->eps_patterned[i], rcwa->G, rcwa->Epsilon2[i], rcwa->Epsilon_inv[i]);
    SolveLayerEigensystem(
			  rcwa->omega_complex,
			  rcwa->nG,
			  rcwa->kx_fwd,
			  rcwa->ky_fwd,
			  rcwa->Epsilon_inv[i],
			  rcwa->Epsilon2[i],
			  rcwa->epstype[rcwa->layerID_patterned[i]],
			  rcwa->q_fwd[rcwa->layerID_patterned[i]],
			  rcwa->kp_fwd[rcwa->layerID_patterned[i]],
			  rcwa->phi_fwd[rcwa->layerID_patterned[i]],
			  NULL,
			  NULL,
			  0
			  );
    SolveLayerEigensystem(
			  rcwa->omega_complex,
			  rcwa->nG,
			  rcwa->kx_adj,
			  rcwa->ky_adj,
			  rcwa->Epsilon_inv[i],
			  rcwa->Epsilon2[i],
			  rcwa->epstype[rcwa->layerID_patterned[i]],
			  rcwa->q_adj[rcwa->layerID_patterned[i]],
			  rcwa->kp_adj[rcwa->layerID_patterned[i]],
			  rcwa->phi_adj[rcwa->layerID_patterned[i]],
			  NULL,
			  NULL,
			  0
			  );
  }
  
  rcwa->qlist_fwd=(const std::complex<double>**)malloc(rcwa->nlayers*sizeof(std::complex<double>*));
  rcwa->kplist_fwd=(const std::complex<double>**)malloc(rcwa->nlayers*sizeof(std::complex<double>*));
  rcwa->philist_fwd=(const std::complex<double>**)malloc(rcwa->nlayers*sizeof(std::complex<double>*));
  rcwa->qlist_adj=(const std::complex<double>**)malloc(rcwa->nlayers*sizeof(std::complex<double>*));
  rcwa->kplist_adj=(const std::complex<double>**)malloc(rcwa->nlayers*sizeof(std::complex<double>*));
  rcwa->philist_adj=(const std::complex<double>**)malloc(rcwa->nlayers*sizeof(std::complex<double>*));
  for(int i=0;i<rcwa->nlayers;i++){
    rcwa->qlist_fwd[i]=rcwa->q_fwd[i];
    rcwa->kplist_fwd[i]=rcwa->kp_fwd[i];
    rcwa->philist_fwd[i]=rcwa->phi_fwd[i];
    rcwa->qlist_adj[i]=rcwa->q_adj[i];
    rcwa->kplist_adj[i]=rcwa->kp_adj[i];
    rcwa->philist_adj[i]=rcwa->phi_adj[i];
  }
  
}

void rcwa_solve_layer(rcwa_ *rcwa, int which_layer){

  //forward fields
  SolveInterior(
		rcwa->nlayers,
		which_layer,
		rcwa->nG,
		rcwa->kx_fwd,rcwa->ky_fwd,
		rcwa->omega_complex,
		rcwa->thickness,
		rcwa->qlist_fwd,
		rcwa->eps_inv_list,
		rcwa->epstype,
		rcwa->kplist_fwd,
		rcwa->philist_fwd,
		rcwa->a0_fwd,
		rcwa->bN_fwd,
		rcwa->ab_internal_fwd,
		NULL,
		NULL,
		0
		);

  //adjoint fields
  SolveInterior(
		rcwa->nlayers,
		which_layer,
		rcwa->nG,
		rcwa->kx_adj,rcwa->ky_adj,
		rcwa->omega_complex,
		rcwa->thickness,
		rcwa->qlist_adj,
		rcwa->eps_inv_list,
		rcwa->epstype,
		rcwa->kplist_adj,
		rcwa->philist_adj,
		rcwa->a0_adj,
		rcwa->bN_adj,
		rcwa->ab_internal_adj,
		NULL,
		NULL,
		0
		);
  
}

void rcwa_translate_and_getfields(rcwa_ *rcwa, int which_layer, double dz_a, double dz_b,  double kbloch[2]){

  TranslateAmplitudes2(
		       rcwa->nG,
		       rcwa->q_fwd[which_layer],
		       dz_a,
		       dz_b,
		       rcwa->ab_internal_fwd
		       );
  
  TranslateAmplitudes2(
		       rcwa->nG,
		       rcwa->q_adj[which_layer],
		       dz_a,
		       dz_b,
		       rcwa->ab_internal_adj
		       );

  size_t ngr[2]={size_t(rcwa->ngrid[0]),size_t(rcwa->ngrid[1])};
  GetFieldOnGrid(
		 rcwa->nG,
		 rcwa->G,
		 rcwa->kx_fwd,rcwa->ky_fwd,
		 rcwa->omega_complex,
		 rcwa->q_fwd[which_layer],
		 rcwa->kp_fwd[which_layer],
		 rcwa->phi_fwd[which_layer],
		 rcwa->eps_inv_list[which_layer],
		 rcwa->epstype[which_layer],
		 rcwa->ab_internal_fwd,
		 ngr,
		 NULL,
		 rcwa->efield_fwd,
		 rcwa->hfield_fwd
		 );  

  GetFieldOnGrid(
		 rcwa->nG,
		 rcwa->G,
		 rcwa->kx_adj,rcwa->ky_adj,
		 rcwa->omega_complex,
		 rcwa->q_adj[which_layer],
		 rcwa->kp_adj[which_layer],
		 rcwa->phi_adj[which_layer],
		 rcwa->eps_inv_list[which_layer],
		 rcwa->epstype[which_layer],
		 rcwa->ab_internal_adj,
		 ngr,
		 NULL,
		 rcwa->efield_adj,
		 rcwa->hfield_adj
		 );

  for(int iy=0;iy<rcwa->ngrid[1];iy++){
    for(int ix=0;ix<rcwa->ngrid[0];ix++){
      int i=ix+rcwa->ngrid[0]*iy;
      double kx=kbloch[0], ky=kbloch[1];
      double x=(ix+0.5)*rcwa->dg[0], y=(iy+0.5)*rcwa->dg[1];
      std::complex<double> blochfac_fwd=std::exp(std::complex<double>(0, kx*x + ky*y));
      std::complex<double> blochfac_adj=std::exp(std::complex<double>(0,-kx*x - ky*y));
      rcwa->efield_fwd[3*i+0] *= blochfac_fwd;
      rcwa->hfield_fwd[3*i+0] *= blochfac_fwd;
      rcwa->efield_fwd[3*i+1] *= blochfac_fwd;
      rcwa->hfield_fwd[3*i+1] *= blochfac_fwd;
      rcwa->efield_fwd[3*i+2] *= blochfac_fwd;
      rcwa->hfield_fwd[3*i+2] *= blochfac_fwd;

      rcwa->efield_adj[3*i+0] *= blochfac_adj;
      rcwa->hfield_adj[3*i+0] *= blochfac_adj;
      rcwa->efield_adj[3*i+1] *= blochfac_adj;
      rcwa->hfield_adj[3*i+1] *= blochfac_adj;
      rcwa->efield_adj[3*i+2] *= blochfac_adj;
      rcwa->hfield_adj[3*i+2] *= blochfac_adj;

      rcwa->efield_fwd[3*i+0] = std::conj(rcwa->efield_fwd[3*i+0]);
      rcwa->hfield_fwd[3*i+0] = std::conj(rcwa->hfield_fwd[3*i+0]);
      rcwa->efield_fwd[3*i+1] = std::conj(rcwa->efield_fwd[3*i+1]);
      rcwa->hfield_fwd[3*i+1] = std::conj(rcwa->hfield_fwd[3*i+1]);
      rcwa->efield_fwd[3*i+2] = std::conj(rcwa->efield_fwd[3*i+2]);
      rcwa->hfield_fwd[3*i+2] = std::conj(rcwa->hfield_fwd[3*i+2]);
      rcwa->efield_adj[3*i+0] = std::conj(rcwa->efield_adj[3*i+0]);
      rcwa->hfield_adj[3*i+0] = std::conj(rcwa->hfield_adj[3*i+0]);
      rcwa->efield_adj[3*i+1] = std::conj(rcwa->efield_adj[3*i+1]);
      rcwa->hfield_adj[3*i+1] = std::conj(rcwa->hfield_adj[3*i+1]);
      rcwa->efield_adj[3*i+2] = std::conj(rcwa->efield_adj[3*i+2]);
      rcwa->hfield_adj[3*i+2] = std::conj(rcwa->hfield_adj[3*i+2]);      

    }
  }
  

}

void rcwa_getfields_fwd(rcwa_ *rcwa, int which_layer, double dz, double kbloch[2]){

  //forward fields

  SolveInterior(
		rcwa->nlayers,
		which_layer,
		rcwa->nG,
		rcwa->kx_fwd,rcwa->ky_fwd,
		rcwa->omega_complex,
		rcwa->thickness,
		rcwa->qlist_fwd,
		rcwa->eps_inv_list,
		rcwa->epstype,
		rcwa->kplist_fwd,
		rcwa->philist_fwd,
		rcwa->a0_fwd,
		rcwa->bN_fwd,
		rcwa->ab_internal_fwd,
		NULL,
		NULL,
		0
		);

  TranslateAmplitudes(
		      rcwa->nG,
		      rcwa->q_fwd[which_layer],
		      rcwa->thickness[which_layer],
		      dz,
		      rcwa->ab_internal_fwd
		      );
  size_t ngr[2]={size_t(rcwa->ngrid[0]),size_t(rcwa->ngrid[1])};
  GetFieldOnGrid(
		 rcwa->nG,
		 rcwa->G,
		 rcwa->kx_fwd,rcwa->ky_fwd,
		 rcwa->omega_complex,
		 rcwa->q_fwd[which_layer],
		 rcwa->kp_fwd[which_layer],
		 rcwa->phi_fwd[which_layer],
		 rcwa->eps_inv_list[which_layer],
		 rcwa->epstype[which_layer],
		 rcwa->ab_internal_fwd,
		 ngr,
		 NULL,
		 rcwa->efield_fwd,
		 rcwa->hfield_fwd
		 );

  for(int iy=0;iy<rcwa->ngrid[1];iy++){
    for(int ix=0;ix<rcwa->ngrid[0];ix++){
      int i=ix+rcwa->ngrid[0]*iy;
      double kx=kbloch[0], ky=kbloch[1];
      double x=(ix+0.5)*rcwa->dg[0], y=(iy+0.5)*rcwa->dg[1];
      std::complex<double> blochfac_fwd=std::exp(std::complex<double>(0, kx*x + ky*y));
      rcwa->efield_fwd[3*i+0] *= blochfac_fwd;
      rcwa->hfield_fwd[3*i+0] *= blochfac_fwd;
      rcwa->efield_fwd[3*i+1] *= blochfac_fwd;
      rcwa->hfield_fwd[3*i+1] *= blochfac_fwd;
      rcwa->efield_fwd[3*i+2] *= blochfac_fwd;
      rcwa->hfield_fwd[3*i+2] *= blochfac_fwd;

      rcwa->efield_fwd[3*i+0] = std::conj(rcwa->efield_fwd[3*i+0]);
      rcwa->hfield_fwd[3*i+0] = std::conj(rcwa->hfield_fwd[3*i+0]);
      rcwa->efield_fwd[3*i+1] = std::conj(rcwa->efield_fwd[3*i+1]);
      rcwa->hfield_fwd[3*i+1] = std::conj(rcwa->hfield_fwd[3*i+1]);
      rcwa->efield_fwd[3*i+2] = std::conj(rcwa->efield_fwd[3*i+2]);
      rcwa->hfield_fwd[3*i+2] = std::conj(rcwa->hfield_fwd[3*i+2]);

    }
  }

}

void destroy_rcwa(rcwa_ *rcwa){

  free(rcwa->layerID_uniform);
  free(rcwa->layerID_patterned);
  free(rcwa->thickness);
  my_free(rcwa->G);
  free(rcwa->eps_uniform);
  free(rcwa->eps_uniform_inv);
  for(int i=0;i<rcwa->nlayers_patterned;i++){
    free(rcwa->eps_patterned[i]);
    my_free(rcwa->Epsilon2[i]);
    my_free(rcwa->Epsilon_inv[i]);
  }
  free(rcwa->eps_patterned);
  free(rcwa->Epsilon2);
  free(rcwa->Epsilon_inv);
  free(rcwa->eps_inv_list);
  free(rcwa->epstype);
  free(rcwa->kx_fwd);
  free(rcwa->ky_fwd);
  free(rcwa->kx_adj);
  free(rcwa->ky_adj);
  for(int i=0;i<rcwa->nlayers;i++){
    my_free(rcwa->q_fwd[i]);
    my_free(rcwa->kp_fwd[i]);
    my_free(rcwa->phi_fwd[i]);
    my_free(rcwa->q_adj[i]);
    my_free(rcwa->kp_adj[i]);
    my_free(rcwa->phi_adj[i]);
  }
  free(rcwa->qlist_fwd);
  free(rcwa->kplist_fwd);
  free(rcwa->philist_fwd);
  free(rcwa->qlist_adj);
  free(rcwa->kplist_adj);
  free(rcwa->philist_adj);
  free(rcwa->q_fwd);
  free(rcwa->kp_fwd);
  free(rcwa->phi_fwd);
  free(rcwa->q_adj);
  free(rcwa->kp_adj);
  free(rcwa->phi_adj);
  my_free(rcwa->ab_fwd);
  my_free(rcwa->ab_adj);
  if(rcwa->ab_internal_fwd) my_free(rcwa->ab_internal_fwd);
  if(rcwa->ab_internal_adj) my_free(rcwa->ab_internal_adj);
  my_free(rcwa->efield_fwd);
  my_free(rcwa->hfield_fwd);
  my_free(rcwa->efield_adj);
  my_free(rcwa->hfield_adj);

}

void set_cell(cell_ *cell,
	      int *layerID_uniform, double *n_uniform,
	      int *layerID_patterned, double *epsdiff_patterned, double *epsbkg_patterned, int *nz_integrate_patterned,
	      double *thickness){

  cell->nlayers=cell->nlayers_uniform+cell->nlayers_patterned;
  cell->layerID_uniform        = (int *)   malloc(cell->nlayers_uniform  *sizeof(int));
  cell->layerID_patterned      = (int *)   malloc(cell->nlayers_patterned*sizeof(int));
  cell->n_uniform              = (double *)malloc(cell->nlayers_uniform  *sizeof(double));
  cell->epsdiff_patterned      = (double *)malloc(cell->nlayers_patterned*sizeof(double));
  cell->epsbkg_patterned       = (double *)malloc(cell->nlayers_patterned*sizeof(double));
  cell->nz_integrate_patterned = (int *)   malloc(cell->nlayers_patterned*sizeof(int));
  cell->thickness              = (double *)malloc(cell->nlayers          *sizeof(double));

  for(int i=0;i<cell->nlayers_uniform;i++){
    cell->layerID_uniform[i]=layerID_uniform[i];
    cell->n_uniform[i]=n_uniform[i];
  }
  for(int i=0;i<cell->nlayers_patterned;i++){
    cell->layerID_patterned[i]=layerID_patterned[i];
    cell->epsdiff_patterned[i]=epsdiff_patterned[i];
    cell->epsbkg_patterned[i]=epsbkg_patterned[i];
    cell->nz_integrate_patterned[i]=nz_integrate_patterned[i];
  }
  for(int i=0;i<cell->nlayers;i++)
    cell->thickness[i]=thickness[i];
  
  int ng2=cell->ngrid[0]*cell->ngrid[1];
  cell->bx=(std::complex<double> *)malloc(ng2*sizeof(std::complex<double>));
  cell->by=(std::complex<double> *)malloc(ng2*sizeof(std::complex<double>));
  cell->vx_cx=(std::complex<double> *)malloc(ng2*sizeof(std::complex<double>));
  cell->vy_cx=(std::complex<double> *)malloc(ng2*sizeof(std::complex<double>));
  cell->vz_cx=(std::complex<double> *)malloc(ng2*sizeof(std::complex<double>));
  cell->vx_cy=(std::complex<double> *)malloc(ng2*sizeof(std::complex<double>));
  cell->vy_cy=(std::complex<double> *)malloc(ng2*sizeof(std::complex<double>));
  cell->vz_cy=(std::complex<double> *)malloc(ng2*sizeof(std::complex<double>));

  cell->vxgrad_complex=(std::complex<double>*)malloc(cell->nlayers_patterned*ng2*sizeof(std::complex<double>));
  cell->vygrad_complex=(std::complex<double>*)malloc(cell->nlayers_patterned*ng2*sizeof(std::complex<double>));
  cell->vzgrad_complex=(std::complex<double>*)malloc(cell->nlayers_patterned*ng2*sizeof(std::complex<double>));

}

void free_cell(cell_ *cell){

  free(cell->layerID_uniform);
  free(cell->layerID_patterned);
  free(cell->n_uniform);
  free(cell->epsdiff_patterned);
  free(cell->epsbkg_patterned);
  free(cell->nz_integrate_patterned);
  free(cell->thickness);
  free(cell->bx);
  free(cell->by);
  free(cell->vx_cx);
  free(cell->vy_cx);
  free(cell->vz_cx);
  free(cell->vx_cy);
  free(cell->vy_cy);
  free(cell->vz_cy);
  free(cell->vxgrad_complex);
  free(cell->vygrad_complex);
  free(cell->vzgrad_complex);

}
