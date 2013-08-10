#ifndef _PLATE_1D_OPTIMIZER_
#define _PLATE_1D_OPTIMIZER_ 1

#include <iostream>
#include <fstream>
#include "Eigen/Eigen"
#include "plate_var_types.h"
#include "Solver.h"

using namespace Eigen;
using std::cout;
using std::endl;
using std::ofstream;

class Optimizer
{
public:
	Optimizer( Solver* _solver, N_PRES _weight, N_PRES _J0h, N_PRES _tauh, N_PRES charTime );
	void optimize( N_PRES Jstart, N_PRES tauStart );

private:
	void calc1stOrdOptInfo( long double J0, long double tau, long double* _objVal, Matrix<N_PRES, 2, 1>* _gk );
	PL_NUM calcFuncVal();
	N_PRES calcBettaN();
	N_PRES calcBettaN_();
	int lineSearch();
	void secant2( N_PRES* _a, N_PRES* _b );
	void update( N_PRES a, N_PRES b, N_PRES c, N_PRES* _a, N_PRES* _b );

	N_PRES min( N_PRES a, N_PRES b );
	N_PRES max( N_PRES a, N_PRES b );

	Solver* solver;

	N_PRES weight;

	N_PRES J0h;
	N_PRES tauh;

	N_PRES charTime;

	Matrix<N_PRES, 2, 1> parVect;
	Matrix<N_PRES, 2, 1> parVectPrev;

	N_PRES objVal;
	N_PRES objValPrev;
	
	Matrix<N_PRES, 2, 1> gk1;		//new gradient
	Matrix<N_PRES, 2, 1> gk;		//gradient from prev step

	Matrix<N_PRES, 2, 1> dk1;		//new descent direction
	Matrix<N_PRES, 2, 1> dk;		//descent direction from previous step

	N_PRES alphak;

	N_PRES linSa;
	N_PRES linSb;
	N_PRES linSc;
	N_PRES linSd;

	N_PRES phi0;
	N_PRES phiPr0;
	N_PRES phiAk;
	N_PRES phiA;
	N_PRES phiPrA;
	N_PRES phiB;
	N_PRES phiPrB;
	N_PRES phiC;
	N_PRES phiPrC;
	N_PRES phiD;
	N_PRES phiPrD;

	N_PRES phiPrAA;
	N_PRES phiPrBB;

	N_PRES wolfeDelta;
	N_PRES wolfeSigma;

	N_PRES epsK;
	N_PRES thetaUpdate;
	N_PRES gammaLineS;

	ofstream dmpo;
};

#endif