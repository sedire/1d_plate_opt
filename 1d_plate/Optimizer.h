#ifndef _PLATE_1D_OPTIMIZER_
#define _PLATE_1D_OPTIMIZER_ 1

#include <iostream>
#include <fstream>
#include "Eigen/Eigen"
#include "plate_var_types.h"
#include "hyperDual.h"
#include "Solver.h"
#include "asa_user.h"
#include "HagerOptFuncs.h"

using namespace Eigen;
using std::cout;
using std::endl;
using std::ofstream;

//double calc1stOrdOptInfoCG_DES( double* g, double* x, long n );

template<class PL_NUM>
void optimizeASA_Taus( const Matrix<N_PRES, GRAD_SIZE_FULL, 1>& params )
{
	time_t totOptStart = time( 0 );
	cout << "optimizeASA enter\n";

	const N_PRES threshold( 1.e-6 * W_SCALE * W_SCALE );

	double* x = new double[GRAD_SIZE_FULL];
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		x[i] = params( i );
	}
	double* lo = new double[GRAD_SIZE_FULL];
	double* hi = new double[GRAD_SIZE_FULL];
	for( int i = 0; i < GRAD_SIZE_FULL; i += 3 )
	{
		lo[i] = -1.0;
		lo[i + 1] = 0.002;
		lo[i + 2] = 0.00001;

		hi[i] = 1.0;
		hi[i + 1] = 1000000000.0;
		hi[i + 2] = 1000000000.0;
	}

	asacg_parm cgParm;
    asa_parm asaParm;
	asa_cg_default( &cgParm );
    asa_default( &asaParm );
    cgParm.PrintParms = TRUE;
    cgParm.PrintLevel = 3;
    asaParm.PrintParms = TRUE;
    asaParm.PrintLevel = 3;

	asa_cg( x, lo, hi, GRAD_SIZE_FULL, NULL, &cgParm, &asaParm, threshold, calcValASA_Taus, calcGradASA_Taus, calc1stOrdOptInfoASA_Taus, 0, 0 );

	cout << "\n\n===============\nASA optimization complete. X is:\n";
	ofstream of( "solution.txt" );
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		of << x[i] << endl;
		cout << x[i] << endl;
	}
	of.close();

	cout << " total optimization time : " << time( 0 ) - totOptStart << endl;

	delete[] x;
	delete[] lo;
	delete[] hi;
}

#endif