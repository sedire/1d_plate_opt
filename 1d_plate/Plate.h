#ifndef _PLATE_1D_PLATE_
#define _PLATE_1D_PLATE_ 1

#include "plate_var_types.h"
#include <math.h>
#include <complex>

using std::complex;

class Plate
{
public:
	PL_NUM E1;				//Young's modulus
	PL_NUM E2;				//Young's modulus
	PL_NUM nu21;			//Poisson's ratio	
	PL_NUM rho;				//composite's density

	PL_NUM sigma_x;			//electric conductivity
	PL_NUM sigma_x_mu;

	N_PRES h;				//thickness of the plate
	N_PRES a;				//width of the plate

	Plate();
	~Plate();
	void loadVals( PL_NUM _E1, PL_NUM _E2, PL_NUM _nu21, PL_NUM _rho, N_PRES _h, PL_NUM _sigma_x, N_PRES _a );
};

#endif