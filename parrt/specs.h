#ifndef GUARD_specs_h
#define GUARD_specs_h

#include "petsc.h"
#include "nlopt.h"
#include "filters.h"
#include "output.h"
#include "scatter.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include "input.h"
#include "near2far.h"
#include "wrap.h"

void create_incident_pw(PetscScalar *bx, PetscScalar *by, double freq, double n_medium, double polar, double azimuth, double j[4], double Lx, double Ly, int Ngrid[2]);

void populate_b_and_c(int colour, int ngrid[2], int ng2_global, int ng2_per_comm, int ncells_per_comm,
		      double Lx, double Ly, int ncells_along_x, int ncells_along_y, int Ngrid[2],
		      double freq, double n_medium, double polar, double azimuth, double input_polarization[4],
		      double oxy[2], int symxy[2], double xyzfar[3],
		      cell_ *cell);

#endif
