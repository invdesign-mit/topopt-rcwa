#ifndef _OBJFUNCS_H_
#define _OBJFUNCS_H_

#include <complex>
#include <cstddef>
#include "wrap.h"

void v_dot_E(double *epsbar, cell_ *cell,
	     std::complex<double> *bx, std::complex<double> *by,
	     std::complex<double> *ux, std::complex<double> *uy, double uobj[2], std::complex<double> *ugrad, std::complex<double> *intface_ugrad,
	     std::complex<double> *vx, std::complex<double> *vy, double vobj[2], std::complex<double> *vgrad, std::complex<double> *intface_vgrad,
	     std::complex<double> *wx, std::complex<double> *wy, double wobj[2], std::complex<double> *wgrad, std::complex<double> *intface_wgrad);

#endif
