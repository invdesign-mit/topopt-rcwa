#ifndef GUARD_near2far_h
#define GUARD_near2far_h

#include "petsc.h"
#include "nlopt.h"
#include "filters.h"
#include "output.h"
#include "scatter.h"

void green3d(std::complex<double> *EH, const double *x,
	     double freq, double eps, double mu,
	     const double *x0, int comp, std::complex<double> f0);

void create_near2far(PetscScalar *vx_cx, PetscScalar *vy_cx, PetscScalar *vz_cx,
		     PetscScalar *vx_cy, PetscScalar *vy_cy, PetscScalar *vz_cy,
		     int nxy[2], double Lxy[2], double oxy[2], int symxy[2],
		     double xyzfar[3],
		     double freq, double eps, double mu);

void create_near2far_2d(PetscScalar *vx_cx, PetscScalar *vy_cx, PetscScalar *vz_cx,
			PetscScalar *vx_cy, PetscScalar *vy_cy, PetscScalar *vz_cy,
			int nxy[2], double Lxy[2], double oxy[2], int symxy[2],
			double xyzfar[3],
			double freq, double eps, double mu);

void create_near2far_2d_phase(PetscScalar *vx_cx, PetscScalar *vy_cx, PetscScalar *vz_cx,
			      PetscScalar *vx_cy, PetscScalar *vy_cy, PetscScalar *vz_cy,
			      int nxy[2], double Lxy[2], double oxy[2], int symxy[2],
			      double xyzfar[3],
			      double freq, double eps, double mu);

#endif
