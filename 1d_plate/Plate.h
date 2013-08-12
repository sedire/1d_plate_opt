#ifndef _PLATE_1D_PLATE_
#define _PLATE_1D_PLATE_ 1

#include "plate_var_types.h"
#include <math.h>
#include <complex>

using std::complex;

template<class PL_NUM>
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

template<class PL_NUM>
Plate<PL_NUM>::Plate()
{

}

template<class PL_NUM>
Plate<PL_NUM>::~Plate()
{

}

template<class PL_NUM>
void Plate<PL_NUM>::loadVals( PL_NUM _E1, PL_NUM _E2, PL_NUM _nu21, PL_NUM _rho, N_PRES _h, PL_NUM _sigma_x, N_PRES _a )
{
	E1 = _E1;
	E2 = _E2;
	nu21 = _nu21;
	rho = _rho;
	h = _h;
	a = _a;
	sigma_x = _sigma_x;
	sigma_x_mu = _sigma_x * 0.000001256l;
}

#endif