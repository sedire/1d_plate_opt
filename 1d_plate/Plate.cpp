#include "Plate.h"

Plate::Plate()
{

}

Plate::~Plate()
{

}

void Plate::loadVals( PL_NUM _E1, PL_NUM _E2, PL_NUM _nu21, PL_NUM _rho, PL_NUM _h, PL_NUM _sigma_x, PL_NUM _a )
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