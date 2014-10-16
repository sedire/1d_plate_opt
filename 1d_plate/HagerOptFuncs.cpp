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
	time_t begin = time( 0 );
	N_PRES dt = DELTA_T;
	N_PRES dy = 0.1524 / ( NODES_Y - 1 );

	cout << "try to calc at\n";
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		cout << x[i] << endl;
	}
	cout << " ====\n";

	HPD<N_PRES, GRAD_SIZE_SECOND> J0begin1;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauBeginSin1;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauBeginExp1;

	HPD<N_PRES, GRAD_SIZE_SECOND> J0begin2[SCEN_NUMBER];
	HPD<N_PRES, GRAD_SIZE_SECOND> tauBeginSin2[SCEN_NUMBER];
	HPD<N_PRES, GRAD_SIZE_SECOND> tauBeginExp2[SCEN_NUMBER];

	HPD<N_PRES, GRAD_SIZE_SECOND> B0begin2;

	J0begin1.elems[0] = x[0];
	J0begin1.elems[1] = 1.0l;
	J0begin1.elems[2] = 0.0l;
	J0begin1.elems[3] = 0.0l;
	J0begin1.elems[4] = 0.0l;
	J0begin1.elems[5] = 0.0l;
	J0begin1.elems[6] = 0.0l;

	tauBeginSin1.elems[0] = x[1];
	tauBeginSin1.elems[1] = 0.0l;
	tauBeginSin1.elems[2] = 1.0l;
	tauBeginSin1.elems[3] = 0.0l;
	tauBeginSin1.elems[4] = 0.0l;
	tauBeginSin1.elems[5] = 0.0l;
	tauBeginSin1.elems[6] = 0.0l;

	tauBeginExp1.elems[0] = x[2];
	tauBeginExp1.elems[1] = 0.0l;
	tauBeginExp1.elems[2] = 0.0l;
	tauBeginExp1.elems[3] = 1.0l;
	tauBeginExp1.elems[4] = 0.0l;
	tauBeginExp1.elems[5] = 0.0l;
	tauBeginExp1.elems[6] = 0.0l;

	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		J0begin2[scen].elems[0] = x[( scen + 1 ) * 3];
		J0begin2[scen].elems[1] = 0.0l;
		J0begin2[scen].elems[2] = 0.0l;
		J0begin2[scen].elems[3] = 0.0l;
		J0begin2[scen].elems[4] = 1.0l;
		J0begin2[scen].elems[5] = 0.0l;
		J0begin2[scen].elems[6] = 0.0l;

		tauBeginSin2[scen].elems[0] = x[( scen + 1 ) * 3 + 1];
		tauBeginSin2[scen].elems[1] = 0.0l;
		tauBeginSin2[scen].elems[2] = 0.0l;
		tauBeginSin2[scen].elems[3] = 0.0l;
		tauBeginSin2[scen].elems[4] = 0.0l;
		tauBeginSin2[scen].elems[5] = 1.0l;
		tauBeginSin2[scen].elems[6] = 0.0l;

		tauBeginExp2[scen].elems[0] = x[( scen + 1 ) * 3 + 2];
		tauBeginExp2[scen].elems[1] = 0.0l;
		tauBeginExp2[scen].elems[2] = 0.0l;
		tauBeginExp2[scen].elems[3] = 0.0l;
		tauBeginExp2[scen].elems[4] = 0.0l;
		tauBeginExp2[scen].elems[5] = 0.0l;
		tauBeginExp2[scen].elems[6] = 1.0l;
	}

	B0begin2 = 1.0l;

	Solver<HPD<N_PRES, GRAD_SIZE_SECOND> > solver_second[SCEN_NUMBER];

	cout << "\tcalculating func val\n";

	double charTime = CHAR_TIME;

	HPD<N_PRES, GRAD_SIZE_SECOND> funcVal1[SCEN_NUMBER];
	HPD<N_PRES, GRAD_SIZE_SECOND> funcVal2[SCEN_NUMBER];
	N_PRES mechLoad[SCEN_NUMBER] = { 7500000, 10000000, 20000000 };
	N_PRES mechTaus[SCEN_NUMBER] = { 0.008, 0.01, 0.012 };

#pragma omp parallel for
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		cout << omp_get_thread_num() << endl;
		
		solver_second[scen].setTask( J0begin1, tauBeginSin1, tauBeginExp1, J0begin2[scen], tauBeginSin2[scen], tauBeginExp2[scen], B0begin2, 10000000, 0.01 );
		solver_second[scen].setMechLoad( mechLoad[scen], mechTaus[scen] );

		HPD<N_PRES, GRAD_SIZE_SECOND> sum = 0.0;
		funcVal1[scen] = 0.0l;
		HPD<N_PRES, GRAD_SIZE_SECOND> val1;

		while( solver_second[scen].cur_t <= SWITCH_TIME )
		{
			sum += solver_second[scen].do_step();

			solver_second[scen].increaseTime();

			//solver_second.dump_check_sol( -1 );
		}
		funcVal1[scen] = sum * dt * dy / SWITCH_TIME;

		funcVal2[scen] = 0.0l;
		HPD<N_PRES, GRAD_SIZE_SECOND> val2;
		sum = 0.0;
		while( solver_second[scen].cur_t <= charTime )
		{
			//cout << "\t\t both -- " << solver.cur_t.real() << " params: " << x[0] << " " << x[1] << " " << x[2] << endl;
			sum += solver_second[scen].do_step();

			solver_second[scen].increaseTime();
			//solver_second[scen].cur_t += solver_second[scen].dt;
			//++( solver_second[i].curTimeStep );

			//solver_second.dump_check_sol( -1 );
		}
		funcVal2[scen] = sum * dt * dy / ( charTime - SWITCH_TIME );
	}

	cout << "\tfunc val done\n";
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		for( int j = 0; j <= GRAD_SIZE_FIRST; ++j )
		{
			cout << "\t1st; scen " << scen << " " << funcVal1[scen].elems[j] << endl;
		}
		for( int j = 0; j <= GRAD_SIZE_SECOND; ++j )
		{
			cout << "\t2nd; scen " << scen << " " << funcVal2[scen].elems[j] << endl;
		}
	}
	cout << " -------------\n";

	N_PRES Weight = J_WEIGHT;
	if( g != 0 )
	{
		g[0] = ( funcVal1[0].elems[1] + funcVal1[1].elems[1] + funcVal1[2].elems[1]
				+ funcVal2[0].elems[1] + funcVal2[1].elems[1] + funcVal2[2].elems[1] ) / 3.0
				/*+ Weight * 2.0l * exp( - 2.0l * SWITCH_TIME / x[2] ) * sin( M_PI * SWITCH_TIME / x[1] )
				* ( 3.0l * x[0] * sin( M_PI * SWITCH_TIME / x[1] ) + exp( SWITCH_TIME / x[2] ) * (
				- exp( -SWITCH_TIME / x[11] ) * x[9] * sin( M_PI * SWITCH_TIME / x[10] )
				- exp( -SWITCH_TIME / x[5] ) * x[3] * sin( M_PI * SWITCH_TIME / x[4] )
				- exp( -SWITCH_TIME / x[8] ) * x[6] * sin( M_PI * SWITCH_TIME / x[7] ) ) )*/;

		g[1] = ( funcVal1[0].elems[2] + funcVal1[1].elems[2] + funcVal1[2].elems[2]
				+ funcVal2[0].elems[2] + funcVal2[1].elems[2] + funcVal2[2].elems[2] ) / 3.0
				/*+ Weight / x[1] / x[1] * 2.0l * exp( -SWITCH_TIME * ( 2.0 / x[2] + 1.0l / x[5] + 1.0l / x[8] + 1.0l / x[11] ) ) * x[0] * M_PI * SWITCH_TIME
				* cos( M_PI * SWITCH_TIME / x[1] ) * ( -3.0l * exp( SWITCH_TIME * ( 1.0l / x[5] + 1.0l / x[8] + 1.0l / x[11] ) ) * x[0] * sin( M_PI * SWITCH_TIME / x[1] )
				+ exp( SWITCH_TIME * ( 1.0l / x[2] + 1.0l / x[5] + 1.0l / x[8] ) ) * x[9] * sin( M_PI * SWITCH_TIME / x[10] ) 
				+ exp( SWITCH_TIME * ( 1.0l / x[11] + 1.0l / x[2] ) ) * ( exp( SWITCH_TIME / x[8] ) * x[3] * sin( M_PI * SWITCH_TIME / x[4] )
				+ exp( SWITCH_TIME / x[5] ) * x[6] * sin( M_PI * SWITCH_TIME / x[7] ) ) )*/;

		g[2] = ( funcVal1[0].elems[3] + funcVal1[1].elems[3] + funcVal1[2].elems[3]
				+ funcVal2[0].elems[3] + funcVal2[1].elems[3] + funcVal2[2].elems[3] ) / 3.0
				/*- 2.0l / x[2] / x[2] * exp( -SWITCH_TIME * ( 1.0 / x[11] + 2.0 / x[2] + 1.0 / x[5] + 1.0 / x[8] ) ) * SWITCH_TIME * Weight * x[0] * sin( M_PI * SWITCH_TIME / x[1] )
				* ( -3.0l * exp( SWITCH_TIME * ( 1.0 / x[11] + 1.0 / x[5] + 1.0 / x[8] ) ) * x[0] * sin( M_PI * SWITCH_TIME / x[1] )
				+ exp( SWITCH_TIME * ( 1.0 / x[2] + 1.0 / x[5] + 1.0 / x[8] ) ) * x[9]  * sin( M_PI * SWITCH_TIME / x[10] )
				+ exp( SWITCH_TIME * ( 1.0 / x[11] + 1.0 / x[2] ) ) * ( exp( SWITCH_TIME / x[8] ) * x[3] * sin( M_PI * SWITCH_TIME / x[4] )
				+ exp( SWITCH_TIME / x[5] ) * x[6] * sin( M_PI * SWITCH_TIME / x[7] ) ) )*/;

		g[3] = funcVal2[0].elems[4] / 3.0
				/*+ 2.0l * exp( -2.0 * SWITCH_TIME / x[5] ) * Weight * sin( M_PI * SWITCH_TIME / x[4] )
				* ( -exp( SWITCH_TIME * ( -1.0 / x[2] + 1.0 / x[5] ) ) * x[0] * sin( M_PI * SWITCH_TIME / x[1] ) + x[3] * sin( M_PI * SWITCH_TIME / x[4] ) )*/;
		g[4] = funcVal2[0].elems[5] / 3.0
				/*+ 2.0l * exp( -SWITCH_TIME / x[5] ) * M_PI * SWITCH_TIME * Weight * x[3] * cos( M_PI * SWITCH_TIME / x[4] )
				* ( exp( -SWITCH_TIME / x[2] ) * x[0] * sin( M_PI * SWITCH_TIME / x[1] ) - exp( -SWITCH_TIME / x[5] ) * x[3] * sin( M_PI * SWITCH_TIME / x[4] ) ) / x[4] / x[4]*/;
		g[5] = funcVal2[0].elems[6] / 3.0
				/*-2.0l * exp( -SWITCH_TIME / x[5] ) * SWITCH_TIME * Weight * x[3] * sin( M_PI * SWITCH_TIME / x[4] )
				* ( exp( -SWITCH_TIME / x[2] ) * x[0] * sin( M_PI * SWITCH_TIME / x[1] ) - exp( -SWITCH_TIME / x[5] ) * x[3] * sin( M_PI * SWITCH_TIME / x[4] ) ) / x[5] / x[5]*/;

		g[6] = funcVal2[1].elems[4] / 3.0
				/*+ 2.0l * exp( -2.0 * SWITCH_TIME / x[8] ) * Weight * sin( M_PI * SWITCH_TIME / x[7] )
				* ( -exp( SWITCH_TIME * ( -1.0 / x[2] + 1.0 / x[8] ) ) * x[0] * sin( M_PI * SWITCH_TIME / x[1] ) + x[6] * sin( M_PI * SWITCH_TIME / x[7] ) )*/;
		g[7] = funcVal2[1].elems[5] / 3.0
				/*+ 2.0l * exp( -SWITCH_TIME / x[8] ) * M_PI * SWITCH_TIME * Weight * x[6] * cos( M_PI * SWITCH_TIME / x[7] )
				* ( exp( -SWITCH_TIME / x[2] ) * x[0] * sin( M_PI * SWITCH_TIME / x[1] ) - exp( -SWITCH_TIME / x[8] ) * x[6] * sin( M_PI * SWITCH_TIME / x[7] ) ) / x[7] / x[7]*/;
		g[8] = funcVal2[1].elems[6] / 3.0
				/*-2.0l * exp( -SWITCH_TIME / x[8] ) * SWITCH_TIME * Weight * x[6] * sin( M_PI * SWITCH_TIME / x[7] ) 
				* ( exp( -SWITCH_TIME / x[2] ) * x[0] * sin( M_PI * SWITCH_TIME / x[1] ) - exp( -SWITCH_TIME / x[8] ) * x[6] * sin( M_PI * SWITCH_TIME / x[7] ) ) / x[8] / x[8]*/;

		g[9] = funcVal2[2].elems[4] / 3.0
				/*+ 2.0l * exp( -2.0 * SWITCH_TIME / x[11] ) * Weight * sin( M_PI * SWITCH_TIME / x[10] )
				* ( -exp( SWITCH_TIME * ( 1.0 / x[11] - 1.0 / x[2] ) ) * x[0] * sin( M_PI * SWITCH_TIME / x[1] ) + x[9] * sin( M_PI * SWITCH_TIME / x[10] ) )*/;
		g[10] = funcVal2[2].elems[5] / 3.0
				/*+ 2.0l * exp( -SWITCH_TIME / x[11] ) * M_PI * SWITCH_TIME * Weight * x[9] * cos( M_PI * SWITCH_TIME / x[10] )
				* ( exp( -SWITCH_TIME / x[2] ) * x[0] * sin( M_PI * SWITCH_TIME / x[1] ) - exp( -SWITCH_TIME / x[11] ) * x[9] * sin( M_PI * SWITCH_TIME / x[10] ) ) / x[10] / x[10]*/;
		g[11] = funcVal2[2].elems[6] / 3.0
				/*-2.0l * exp( -SWITCH_TIME / x[11] ) * SWITCH_TIME * Weight * x[9] * sin( M_PI * SWITCH_TIME / x[10] ) 
				* ( exp( -SWITCH_TIME / x[2] ) * x[0] * sin( M_PI * SWITCH_TIME / x[1] ) - exp( -SWITCH_TIME / x[11] ) * x[9] * sin( M_PI * SWITCH_TIME / x[10] ) ) / x[11] / x[11]*/;
	}
	ret = ( funcVal1[0].real() + funcVal1[1].real() + funcVal1[2].real()
			+ funcVal2[0].real() + funcVal2[1].real() + funcVal2[2].real() ) / 3.0l
			/*+ Weight * ( ( x[0] * exp( -SWITCH_TIME / x[2] ) * sin( M_PI * SWITCH_TIME / x[1] ) - x[3] * exp( -SWITCH_TIME / x[5] ) * sin( M_PI * SWITCH_TIME / x[4] ) ) *
						( x[0] * exp( -SWITCH_TIME / x[2] ) * sin( M_PI * SWITCH_TIME / x[1] ) - x[3] * exp( -SWITCH_TIME / x[5] ) * sin( M_PI * SWITCH_TIME / x[4] ) ) + 
						( x[0] * exp( -SWITCH_TIME / x[2] ) * sin( M_PI * SWITCH_TIME / x[1] ) - x[6] * exp( -SWITCH_TIME / x[8] ) * sin( M_PI * SWITCH_TIME / x[7] ) ) *
						( x[0] * exp( -SWITCH_TIME / x[2] ) * sin( M_PI * SWITCH_TIME / x[1] ) - x[6] * exp( -SWITCH_TIME / x[8] ) * sin( M_PI * SWITCH_TIME / x[7] ) ) +
						( x[0] * exp( -SWITCH_TIME / x[2] ) * sin( M_PI * SWITCH_TIME / x[1] ) - x[9] * exp( -SWITCH_TIME / x[11] ) * sin( M_PI * SWITCH_TIME / x[10] ) ) * 
						( x[0] * exp( -SWITCH_TIME / x[2] ) * sin( M_PI * SWITCH_TIME / x[1] ) - x[9] * exp( -SWITCH_TIME / x[11] ) * sin( M_PI * SWITCH_TIME / x[10] ) ) )*/;

	time_t endtime = time( 0 );
	cout << "\tdone in " << endtime - begin << endl;

	return ret;
}

double calcValTaus( double* x, long n )
{
	double ret = 0.0l;
	time_t begin = time( 0 );
	N_PRES dt = DELTA_T;
	N_PRES dy = 0.1524 / ( NODES_Y - 1 );

	cout << "try to calc at\n";
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		cout << x[i] << endl;
	}
	cout << " ====\n";

	N_PRES J0begin;
	N_PRES tauBeginSin;
	N_PRES tauBeginExp;

	N_PRES J0begin2[SCEN_NUMBER];
	N_PRES tauBeginSin2[SCEN_NUMBER];
	N_PRES tauBeginExp2[SCEN_NUMBER];

	N_PRES B0begin;

	J0begin = x[0];
	tauBeginSin = x[1];
	tauBeginExp = x[2];

	for( int i = 0; i < SCEN_NUMBER; ++i )
	{
		J0begin2[i] = x[( i + 1 ) * 3];
		tauBeginSin2[i] = x[( i + 1 ) * 3 + 1];
		tauBeginExp2[i] = x[( i + 1 ) * 3 + 2];
	}
	B0begin = 1.0l;

	Solver<N_PRES> solver[SCEN_NUMBER];

	cout << "\tcalculating func val\n";

	double charTime = CHAR_TIME;

	N_PRES funcVal1[SCEN_NUMBER];
	N_PRES funcVal2[SCEN_NUMBER];
	N_PRES mechLoad[SCEN_NUMBER] = { 7500000, 10000000, 20000000 };
	N_PRES mechTaus[SCEN_NUMBER] = { 0.008, 0.01, 0.012 };

#pragma omp parallel for
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		funcVal1[scen] = 0.0l;
		funcVal2[scen] = 0.0l;
		solver[scen].setTask( J0begin, tauBeginSin, tauBeginExp, J0begin2[scen], tauBeginSin2[scen], tauBeginExp2[scen], B0begin, 10000000, 0.01 );
		solver[scen].setMechLoad( mechLoad[scen], mechTaus[scen] );
		N_PRES sum = 0.0;

		while( solver[scen].cur_t <= SWITCH_TIME )
		{
			//cout << "\t\t both -- " << solver[scen].cur_t << " params: " << x[0] << " " << x[1] << " " << x[2] << endl;
			sum += solver[scen].do_step();

			solver[scen].increaseTime(); 

			//solver_second.dump_check_sol( -1 );
		}
		funcVal1[scen] = sum * dt * dy / SWITCH_TIME;

		cout << " scen2 " << scen << endl;

		sum = 0.0;
		while( solver[scen].cur_t <= charTime )
		{
			sum += solver[scen].do_step();

			solver[scen].increaseTime();

			//solver_second.dump_check_sol( -1 );
		}
		funcVal2[scen] = sum * dt * dy / ( charTime - SWITCH_TIME );
	}

	cout << "\tfunc val done\n";
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		cout << "\t1st; scen " << scen << " " << funcVal1[scen] << endl;
		cout << "\t2nd; scen " << scen << " " << funcVal2[scen] << endl;
	}
	cout << " -------------\n";

	N_PRES Weight = J_WEIGHT;
	ret = ( funcVal1[0] + funcVal1[1] + funcVal1[2]
		+ funcVal2[0] + funcVal2[1] + funcVal2[2] ) / 3.0
			/*+ Weight * ( ( x[0] * exp( -SWITCH_TIME / x[2] ) * sin( M_PI * SWITCH_TIME / x[1] ) - x[3] * exp( -SWITCH_TIME / x[5] ) * sin( M_PI * SWITCH_TIME / x[4] ) ) *
						( x[0] * exp( -SWITCH_TIME / x[2] ) * sin( M_PI * SWITCH_TIME / x[1] ) - x[3] * exp( -SWITCH_TIME / x[5] ) * sin( M_PI * SWITCH_TIME / x[4] ) ) + 
						( x[0] * exp( -SWITCH_TIME / x[2] ) * sin( M_PI * SWITCH_TIME / x[1] ) - x[6] * exp( -SWITCH_TIME / x[8] ) * sin( M_PI * SWITCH_TIME / x[7] ) ) *
						( x[0] * exp( -SWITCH_TIME / x[2] ) * sin( M_PI * SWITCH_TIME / x[1] ) - x[6] * exp( -SWITCH_TIME / x[8] ) * sin( M_PI * SWITCH_TIME / x[7] ) ) +
						( x[0] * exp( -SWITCH_TIME / x[2] ) * sin( M_PI * SWITCH_TIME / x[1] ) - x[9] * exp( -SWITCH_TIME / x[11] ) * sin( M_PI * SWITCH_TIME / x[10] ) ) * 
						( x[0] * exp( -SWITCH_TIME / x[2] ) * sin( M_PI * SWITCH_TIME / x[1] ) - x[9] * exp( -SWITCH_TIME / x[11] ) * sin( M_PI * SWITCH_TIME / x[10] ) ) )*/;

	time_t endtime = time( 0 );
	cout << "\tdone in " << endtime - begin << endl;

	cout << funcVal1[0] + funcVal2[0] << endl;
	cout << funcVal1[1] + funcVal2[1] << endl;
	cout << funcVal1[2] + funcVal2[2] << endl;

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
	cout << "\tcalc 1st order CG_DES Both\n";
	double ret = calcValGradTaus( asa->g, asa->x, asa->n );
	return ret;
}

double calcValASA_Taus( asa_objective* asa )
{
	cout << "\tcalc 1st order CG_DES Val\n";
	//return calcValGradTaus( 0, asa->x, asa->n );
	return calcValTaus( asa->x, asa->n );
}

void calcGradASA_Taus( asa_objective* asa )
{
	cout << "\tcalc 1st order CG_DES Grad\n";
	calcValGradTaus( asa->g, asa->x, asa->n );
}
