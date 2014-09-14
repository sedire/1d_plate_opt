#ifndef _PLATE_1D_SOLVERPAR_
#define _PLATE_1D_SOLVERPAR_ 1

#include "plate_var_types.h"

struct SolverPar
{
//parameters of the material
	N_PRES E1;				//Young's modulus
	N_PRES E2;				//Young's modulus
	N_PRES nu21;			//Poisson's ratio	
	N_PRES rho;				//composite's density

	N_PRES sigma_x;			//electric conductivity
	N_PRES sigma_x_mu;

	N_PRES h;				//thickness of the plate
	N_PRES a;				//width of the plate
//scheme params
	N_PRES dt;
	int Km;
	N_PRES dx;
//stress
	N_PRES J0;
	N_PRES J0_1;
	N_PRES tauSin;
	N_PRES tauSin_1;
	N_PRES tauExp;
	N_PRES tauExp_1;

	N_PRES  p0;				//constant mechanical load
	N_PRES tauP;
	N_PRES rad;

	int stress_type;
	int current_type;
//Newmark params
	N_PRES beta;
//some other
	N_PRES B11;
	N_PRES B22;
	N_PRES B12;
	N_PRES By0;
	N_PRES By1;                                      // in considered boundary-value problem

	N_PRES eps_0;
	N_PRES eps_x;
	N_PRES eps_x_0;
};

#endif