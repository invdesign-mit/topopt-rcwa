#ifndef _WRAP_H_
#define _WRAP_H_

#include <complex>
#include <cstddef>

typedef struct {

  int nlayers_uniform;
  int nlayers_patterned;
  int nlayers;
  int *layerID_uniform;
  int *layerID_patterned;
  double *thickness;
  unsigned int nG;
  int *G;
  double Lr[4];
  double Lk[4];
  int ngrid[2];
  double dg[2];
  std::complex<double> *eps_uniform;
  std::complex<double> *eps_uniform_inv;
  std::complex<double> **eps_patterned;
  double freq;
  double omega;
  std::complex<double> omega_complex;
  double kbloch[2];
  double *kx_fwd;
  double *ky_fwd;
  double *kx_adj;
  double *ky_adj;
  std::complex<double> **Epsilon2;
  std::complex<double> **Epsilon_inv;
  int *epstype;
  const std::complex<double> **eps_inv_list;
  std::complex<double> **q_fwd;
  std::complex<double> **kp_fwd;
  std::complex<double> **phi_fwd;
  std::complex<double> **q_adj;
  std::complex<double> **kp_adj;
  std::complex<double> **phi_adj;
  const std::complex<double> **qlist_fwd;
  const std::complex<double> **kplist_fwd;
  const std::complex<double> **philist_fwd;
  const std::complex<double> **qlist_adj;
  const std::complex<double> **kplist_adj;
  const std::complex<double> **philist_adj;
  std::complex<double> *ab_fwd;
  std::complex<double> *a0_fwd;
  std::complex<double> *bN_fwd;
  std::complex<double> *a0x_fwd;
  std::complex<double> *a0y_fwd;
  std::complex<double> *bNx_fwd;
  std::complex<double> *bNy_fwd;
  std::complex<double> *ab_adj;
  std::complex<double> *a0_adj;
  std::complex<double> *bN_adj;
  std::complex<double> *a0x_adj;
  std::complex<double> *a0y_adj;
  std::complex<double> *bNx_adj;
  std::complex<double> *bNy_adj;
  std::complex<double> *ab_internal_fwd;
  std::complex<double> *ab_internal_adj;
  std::complex<double> *efield_fwd;
  std::complex<double> *hfield_fwd;
  std::complex<double> *efield_adj;
  std::complex<double> *hfield_adj;
  
} rcwa_;

typedef struct {

  int nlayers_uniform;
  int nlayers_patterned;
  int nlayers;
  int *layerID_uniform;
  int *layerID_patterned;
  double *thickness;
  unsigned int nG;
  double Lr[4];
  int ngrid[2];
  double freq;
  double *n_uniform;
  double dir_cosine[2];
  double *epsdiff_patterned;
  double *epsbkg_patterned;
  int *nz_integrate_patterned; // # of discretized z points for PATTERNED layers; choose even numbers for Simpson's rule

  std::complex<double> *bx;
  std::complex<double> *by;
  std::complex<double> *vx_cx;
  std::complex<double> *vy_cx;
  std::complex<double> *vz_cx;
  std::complex<double> *vx_cy;
  std::complex<double> *vy_cy;
  std::complex<double> *vz_cy;
  std::complex<double> *vxgrad_complex;
  std::complex<double> *vygrad_complex;
  std::complex<double> *vzgrad_complex;
  
} cell_;

void memset_rcwa(rcwa_ *rcwa, double *epsbar, double kbloch[2],
		 int *layerID_uniform, double *n_uniform,
		 int *layerID_patterned, double *epsdiff, double *epsbkg,
		 double *thickness);

void rcwa_solve_layer(rcwa_ *rcwa, int which_layer);

void rcwa_translate_and_getfields(rcwa_ *rcwa, int which_layer, double dz_a, double dz_b,  double kbloch[2]);

void rcwa_getfields_fwd(rcwa_ *rcwa, int which_layer, double dz, double kbloch[2]);

void destroy_rcwa(rcwa_ *rcwa);

void set_cell(cell_ *cell,
	      int *layerID_uniform, double *n_uniform,
	      int *layerID_patterned, double *epsdiff_patterned, double *epsbkg_patterned, int *nz_integrate_patterned,
	      double *thickness);

void free_cell(cell_ *cell);

#endif
