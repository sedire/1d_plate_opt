#include "plate_var_types.h"
#include "HagerOptFuncs.h"
#include "hyperDual.h"
#include "time.h"
#include "Plate.h"
#include "Solver.h"
#include <iostream>

using std::cout;
using std::endl;

//double calc1stOrdOptInfoCG_DES( double* g, double* x, long n )
//{
//	double ret = 0.0l;
//
//	cout << "\tcalc 1st order CG_DES Both\n";
//	time_t begin = time( 0 );
//
//	HPD<N_PRES, GRAD_SIZE> J0begin;
//	HPD<N_PRES, GRAD_SIZE> tauBegin;
//	HPD<N_PRES, GRAD_SIZE> tauBeginSin;
//	HPD<N_PRES, GRAD_SIZE> B0begin;
//
//	J0begin.elems[0] = x[0];
//	J0begin.elems[1] = 1.0l;
//	J0begin.elems[2] = 0.0l;
//	J0begin.elems[3] = 0.0l;
//
//	tauBegin.elems[0] = x[1];
//	tauBegin.elems[1] = 0.0l;
//	tauBegin.elems[2] = 1.0l;
//	tauBegin.elems[3] = 0.0l;
//
//	tauBeginSin.elems[0] = x[1];
//	tauBeginSin.elems[1] = 0.0l;
//	tauBeginSin.elems[2] = 1.0l;
//	tauBeginSin.elems[3] = 0.0l;
//
//	B0begin.elems[0] = x[2];
//	B0begin.elems[1] = 0.0l;
//	B0begin.elems[2] = 0.0l;
//	B0begin.elems[3] = 1.0l;
//
//	Solver<HPD<N_PRES, GRAD_SIZE> > solver;
//	solver.setTask( J0begin, tauBegin, tauBeginSin, B0begin, 10000000, 0.01 );
//
//	cout << "\tcalculating func val\n";
//
//	double charTime = CHAR_TIME;
//	N_PRES weightJ = 0.0;//1000.0;
//	N_PRES weightB = 0;//1.0l / 6.0 / 6.0 / 6.0;
//	HPD<N_PRES, GRAD_SIZE> funcVal;
//
//	while( solver.cur_t <= charTime )
//	{
//		//cout << "\t\t both -- " << solver.cur_t.real() << " params: " << x[0] << " " << x[1] << " " << x[2] << endl;
//		HPD<N_PRES, GRAD_SIZE> val;
//		val = solver.do_step();
//		funcVal += val * val;
//
//		solver.cur_t += solver.dt;
//		++( solver.curTimeStep );
//
//		solver.dump_check_sol( -1 );
//	}
//
//	cout << "\tfunc val done " << funcVal << endl;
//
//	//( *_gk )( 0 ) = funcVal.elems[1] + weightJ * curVal( 0 ) / sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) );
//	//( *_gk )( 1 ) = funcVal.elems[2];
//	//( *_gk )( 2 ) = funcVal.elems[3] + weightJ * curVal( 2 ) / sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) );
//	//*_objVal = funcVal.real() + sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) ) * weightJ;
//
//	if( g != 0 )
//	{
//		g[0] = funcVal.elems[1] + weightJ * x[0] / 2.0l / ( 1.0l + M_PI * M_PI )
//					* x[1] * exp( -2.0l * charTime / x[1] ) * ( M_PI * M_PI * ( exp( 2.0l * charTime / x[1] ) - 1 )
//					- M_PI * sin( 2.0l * M_PI * charTime / x[1] ) + cos( 2.0l * M_PI * charTime / x[1] ) - 1.0l );
//		g[1] = funcVal.elems[2] + weightJ * x[0] * x[0] / 4.0l  / ( 1.0l + M_PI * M_PI ) / x[1]
//					* exp( -2.0l * charTime / x[1] ) * ( M_PI * M_PI * x[1] * exp( 2.0l * charTime / x[1] )
//					- ( 1.0l + M_PI * M_PI ) * ( x[1] + 2.0l * charTime ) - M_PI * x[1] * sin( 2.0l * M_PI * charTime / x[1] )
//					+ ( x[1] + 2.0l * ( 1.0l + M_PI * M_PI ) * charTime ) * cos( 2.0l * M_PI * charTime / x[1] ) );
//		g[2] = funcVal.elems[3] + 2.0 * weightB * x[2];
//	}
//	ret = funcVal.real() + weightB * x[2] * x[2] + weightJ * x[0] * x[0] / 4.0l / ( 1.0l + M_PI * M_PI )
//			* x[1] * exp( -2.0l * charTime / x[1] ) * ( M_PI * M_PI * ( exp( 2.0l * charTime / x[1] ) - 1 )
//			- M_PI * sin( 2.0l * M_PI * charTime / x[1] ) + cos( 2.0l * M_PI * charTime / x[1] ) - 1.0l );
//
//	time_t endtime = time( 0 );
//	cout << "\tdone in " << endtime - begin << endl;
//
//	return ret;
//}

double calcValGradTaus( double* g, double* x, long n )
{
	double ret = 0.0l;

	cout << "\tcalc 1st order CG_DES Both\n";
	time_t begin = time( 0 );

	HPD<N_PRES, GRAD_SIZE> J0begin;
	HPD<N_PRES, GRAD_SIZE> tauBeginSin;
	HPD<N_PRES, GRAD_SIZE> tauBeginExp;
	HPD<N_PRES, GRAD_SIZE> J0begin_1;
	HPD<N_PRES, GRAD_SIZE> tauBeginSin_1;
	HPD<N_PRES, GRAD_SIZE> tauBeginExp_1;
	HPD<N_PRES, GRAD_SIZE> B0begin;

	J0begin.elems[0] = x[0];
	J0begin.elems[1] = 1.0l;
	J0begin.elems[2] = 0.0l;
	J0begin.elems[3] = 0.0l;
	J0begin.elems[4] = 0.0l;
	J0begin.elems[5] = 0.0l;
	J0begin.elems[6] = 0.0l;

	tauBeginSin.elems[0] = x[1];
	tauBeginSin.elems[1] = 0.0l;
	tauBeginSin.elems[2] = 1.0l;
	tauBeginSin.elems[3] = 0.0l;
	tauBeginSin.elems[4] = 0.0l;
	tauBeginSin.elems[5] = 0.0l;
	tauBeginSin.elems[6] = 0.0l;

	tauBeginExp.elems[0] = x[2];
	tauBeginExp.elems[1] = 0.0l;
	tauBeginExp.elems[2] = 0.0l;
	tauBeginExp.elems[3] = 1.0l;
	tauBeginExp.elems[4] = 0.0l;
	tauBeginExp.elems[5] = 0.0l;
	tauBeginExp.elems[6] = 0.0l;

	J0begin_1.elems[0] = x[3];
	J0begin_1.elems[1] = 0.0l;
	J0begin_1.elems[2] = 0.0l;
	J0begin_1.elems[3] = 0.0l;
	J0begin_1.elems[4] = 1.0l;
	J0begin_1.elems[5] = 0.0l;
	J0begin_1.elems[6] = 0.0l;

	tauBeginSin_1.elems[0] = x[4];
	tauBeginSin_1.elems[1] = 0.0l;
	tauBeginSin_1.elems[2] = 0.0l;
	tauBeginSin_1.elems[3] = 0.0l;
	tauBeginSin_1.elems[4] = 0.0l;
	tauBeginSin_1.elems[5] = 1.0l;
	tauBeginSin_1.elems[6] = 0.0l;

	tauBeginExp_1.elems[0] = x[5];
	tauBeginExp_1.elems[1] = 0.0l;
	tauBeginExp_1.elems[2] = 0.0l;
	tauBeginExp_1.elems[3] = 0.0l;
	tauBeginExp_1.elems[4] = 0.0l;
	tauBeginExp_1.elems[5] = 0.0l;
	tauBeginExp_1.elems[6] = 1.0l;

	B0begin.elems[0] = 1.0;
	B0begin.elems[1] = 0.0l;
	B0begin.elems[2] = 0.0l;
	B0begin.elems[3] = 0.0l;
	B0begin.elems[4] = 0.0l;
	B0begin.elems[5] = 0.0l;
	B0begin.elems[6] = 0.0l;

	Solver<HPD<N_PRES, GRAD_SIZE> > solver;

	cout << "\tcalculating func val\n";

	double charTime = CHAR_TIME;
	N_PRES weightJ = 0.0;//1000.0;

	HPD<N_PRES, GRAD_SIZE> funcVal;
	N_PRES mechLoad[1] = { /*5000000,*/ 10000000/*, 20000000*/ };
	N_PRES mechTaus[1] = { /*0.008,*/ 0.01/*, 0.012*/ };


	solver.setTask( J0begin, tauBeginSin, tauBeginExp, J0begin_1, tauBeginSin_1, tauBeginExp_1, B0begin, 10000000, 0.01 );
	//solver.setMechLoad( mechLoad[i], mechTaus[i] );
	while( solver.cur_t <= charTime )
	{
		//cout << "\t\t both -- " << solver.cur_t.real() << " params: " << x[0] << " " << x[1] << " " << x[2] << endl;
		HPD<N_PRES, GRAD_SIZE> val;
		val = solver.do_step();
		funcVal += val * val;

		solver.cur_t += solver.dt;
		++( solver.curTimeStep );

		solver.dump_check_sol( -1 );
	}

	cout << "\tfunc val done\n";
	for( int i = 0; i <= GRAD_SIZE; ++i )
	{
		cout << "\t" << funcVal.elems[i] << endl;
	}
	cout << " -------------\n";

	if( g != 0 )
	{
		g[0] = funcVal.elems[1]
			+ weightJ * exp( -2.0 * charTime / x[2] ) * x[0] * x[2] / ( 2.0 * ( M_PI * M_PI * x[2] * x[2] + x[1] * x[1] ) )
			* ( ( -1.0 + exp( 2 * charTime / x[2] ) ) * M_PI * M_PI * x[2] * x[2]  - x[1] * x[1] + x[1] * x[1] * cos( 2.0 * M_PI * charTime / x[1] )
			- M_PI * x[2] * x[1] * sin( 2.0 * M_PI * charTime / x[1] ));
		g[1] = funcVal.elems[2]
			+ weightJ * ( exp( -2.0 * charTime / x[2] ) * x[0] * x[0] * M_PI  * x[2] * ( 2.0 * M_PI * x[2] * ( M_PI * M_PI * charTime * x[2] * x[2]
			+ ( charTime + x[2] ) * x[1] * x[1] ) * cos( 2.0 * M_PI * charTime / x[1] ) + x[1] * ( -2.0 * exp( 2.0 * charTime / x[2] ) * M_PI * x[2] * x[2] * x[1]
			+ ( M_PI * M_PI  * ( 2.0 * charTime - x[2] ) * x[2] * x[2] + ( 2.0 * charTime + x[2] ) * x[1] * x[1] ) * sin( 2.0 * M_PI * charTime / x[1] ) ) ) )
			/ ( 4.0 * x[1] * ( M_PI * M_PI * x[2] * x[2] + x[1] * x[1] ) * ( M_PI * M_PI * x[2] * x[2] + x[1] * x[1] ) );
		g[2] = funcVal.elems[3]
			+ weightJ * ( exp( -2.0 * charTime / x[2] ) * x[0] * x[0] * ( -( 2.0 * charTime + x[2] ) * ( M_PI * M_PI * x[2] * x[2] + x[1] * x[1] )
			* ( M_PI * M_PI * x[2] * x[2] + x[1] * x[1] ) + exp( 2.0 * charTime / x[2] ) * M_PI * M_PI * x[2] * x[2] * x[2] * ( M_PI * M_PI * x[2] * x[2] + 3.0 * x[1] * x[1] )
			+ x[1] * x[1] * ( M_PI * M_PI * ( 2.0 * charTime - x[2] ) * x[2] * x[2] + ( 2.0 * charTime + x[2] ) * x[1] * x[1] ) * cos( 2.0 * M_PI * charTime / x[1] )
			- 2.0 * M_PI * x[2] * x[1] * ( M_PI * M_PI * charTime * x[2] * x[2] + ( charTime + x[2] ) * x[1] * x[1] ) * sin( 2.0 * M_PI * charTime / x[1] ) ) )
			/ ( 4.0 * x[2] * ( M_PI * M_PI * x[2] * x[2] + x[1] * x[1] ) * ( M_PI * M_PI * x[2] * x[2] + x[1] * x[1] ) );
		g[3] = funcVal.elems[4];
		g[4] = funcVal.elems[5];
		g[5] = funcVal.elems[6];
	}
	ret = funcVal.real()
		+ weightJ * exp( -2.0 * charTime / x[2] ) * x[0] * x[0] * x[2] / ( 4.0 * ( M_PI * M_PI * x[2] * x[2] + x[1] * x[1] ) )
		* ( ( -1.0 + exp( 2 * charTime / x[2] ) ) * M_PI * M_PI * x[2] * x[2]  - x[1] * x[1] + x[1] * x[1] * cos( 2.0 * M_PI * charTime / x[1] )
		- M_PI * x[2] * x[1] * sin( 2.0 * M_PI * charTime / x[1] ));

	time_t endtime = time( 0 );
	cout << "\tdone in " << endtime - begin << endl;

	return ret;
}

//double calc1stOrdOptInfoASA( asa_objective* asa )
//{
//	double ret = calc1stOrdOptInfoCG_DES( asa->g, asa->x, asa->n );
//	return ret;
//}
//
//double calcValASA( asa_objective* asa )
//{
//	cout << "\tcalc 1st order CG_DES Val\n";
//	return calc1stOrdOptInfoCG_DES( 0, asa->x, asa->n );
//}
//
//void calcGradASA( asa_objective* asa )
//{
//	cout << "\tcalc 1st order CG_DES Grad\n";
//	calc1stOrdOptInfoCG_DES( asa->g, asa->x, asa->n );
//}

///////////////

double calc1stOrdOptInfoASA_Taus( asa_objective* asa )
{
	double ret = calcValGradTaus( asa->g, asa->x, asa->n );
	return ret;
}

double calcValASA_Taus( asa_objective* asa )
{
	cout << "\tcalc 1st order CG_DES Val\n";
	return calcValGradTaus( 0, asa->x, asa->n );
}

void calcGradASA_Taus( asa_objective* asa )
{
	cout << "\tcalc 1st order CG_DES Grad\n";
	calcValGradTaus( asa->g, asa->x, asa->n );
}
