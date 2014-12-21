#include <iostream>
#include <fstream>
#include <stdio.h>
#include "Solver.h"
#include "Plate.h"
#include "Optimizer.h"
#include "HagerOptFuncs.h"
#include <omp.h>
#include "SolverPar.h"
#include "AdjSolver.h"

using std::cout;
using std::endl;
using std::ifstream;

N_PRES* GlobalResArrays = new N_PRES[SCEN_NUMBER * ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM];
N_PRES* GlobalResDtArrays = new N_PRES[SCEN_NUMBER * ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM];
N_PRES* GlobalResAdjArrays = new N_PRES[SCEN_NUMBER * ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM];

int main()
{
	cout << "hello\n";

	omp_set_num_threads( THREAD_NUM );

	time_t adjTime = 0;
	time_t hpdTime = 0;

	/*N_PRES J0start =  0.01;
	N_PRES tauStart = 0.0048;
	N_PRES tauStartExp = 0.0048;

	N_PRES J0start_1 = 0.01;
	N_PRES tauStart_1 = 0.0048;
	N_PRES tauStartExp_1 = 0.0048;

	N_PRES J0start_2 =  0.01;
	N_PRES tauStart_2 = 0.0048;
	N_PRES tauStartExp_2 = 0.0048;

	N_PRES J0start_3 =  0.01;
	N_PRES tauStart_3 = 0.0048;
	N_PRES tauStartExp_3 = 0.0048;*/

//	N_PRES J0start =  0.0168894;
//	N_PRES tauStart = 0.0100969;
//	N_PRES tauStartExp = 0.134501;
//
//	N_PRES J0start_1 = 1;
//	N_PRES tauStart_1 = 0.00515453;
//	N_PRES tauStartExp_1 = 0.00267568;
//
//	N_PRES J0start_2 =  0.015036;
//	N_PRES tauStart_2 = 0.00586419;
//	N_PRES tauStartExp_2 = 0.00222149;
//
//	N_PRES J0start_3 =  1;
//	N_PRES tauStart_3 = 0.00371065;
//	N_PRES tauStartExp_3 = 0.0030765;
//
//	N_PRES ByStart = 1.0;
//
//////////////////////////////////////
//	Solver</*HPD<*/N_PRES/*, GRAD_SIZE>*/ >* solver = new Solver</*HPD<*/N_PRES/*, GRAD_SIZE>*/ >();
//	solver->setTask( J0start, tauStart, tauStartExp, J0start_2, tauStart_2, tauStartExp_2, ByStart, stress_centered, GlobalP02, GlobalTauP2 );
//	time_t tBegin = time( 0 );
//
//	while( solver->cur_t <= CHAR_TIME )
//	{
//		cout << solver->cur_t << endl;
//		solver->do_step();
//		solver->dumpCheckSol( -1, 0 );
//		//solver->dumpWholeSol( 1 );
//
//		solver->increaseTime();
//	}
//	if( solver->getMaxNewtonIterReached() == 1 )
//	{
//		cout << "max newton iter reached\n";
//	}
//
//	time_t tEnd = time( 0 );
//	cout << " \n computations are done in " << tEnd - tBegin << endl;
//	cout << ".........\n";
//	cout << "... done!\n";
//
//
//	cout << " adjTime " << adjTime << endl;
//	cout << " hpdTIme " << hpdTime << endl;

/////////////////////////////////////////

	//Matrix<N_PRES, GRAD_SIZE_FULL, 1> params;
	//params << J0start, tauStart, tauStartExp, 
	//	J0start_1, tauStart_1, tauStartExp_1,
	//	J0start_2, tauStart_2, tauStartExp_2,
	//	J0start_3, tauStart_3, tauStartExp_3;

	//double x[GRAD_SIZE_FULL] = { 0.0267403, 0.0110724, 0.0213401, 
	//							1.0, 0.00514927, 0.00271984,
	//							0.00991937, 0.00682416, 0.00001,
	//							1.0, 0.00371071, 0.00313839 };
	double x[GRAD_SIZE_FULL] = { 0.01, 0.0048, 0.0048, 
								0.01, 0.0048, 0.0048, 
								0.01, 0.0048, 0.0048, 
								0.01, 0.0048, 0.0048 };

	double gAdj[GRAD_SIZE_FULL];
	double g[GRAD_SIZE_FULL];

	double valAdj = calcValGradTausAdj( gAdj, x, 0 );
	cout << " adj comput done\n";
	double val = 0.0;//calcValGradTaus( g, x, 0 );

	ifstream iff( "HPDgrad.txt" );
	if( iff.is_open() )
	{
		for( int i = 0; i < GRAD_SIZE_FULL; ++i )
		{
			iff >> g[i];
		}
		iff.close();
	}

	cout << " ----------\n";
	cout << " == " << val << " " << valAdj << endl;
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		cout << " :: " << g[i] << " " << gAdj[i] << " " << fabs( ( g[i] - gAdj[i] ) / g[i] ) * 100.0 << " % " << endl;
	}
	cout << " ----------\n";

	//optimizeASA_Taus<HPD<N_PRES, GRAD_SIZE> >( params );

	cout << "\n -- Deleting the solution arrays now...\n";

	delete[] GlobalResArrays;
	delete[] GlobalResDtArrays;
	delete[] GlobalResAdjArrays;

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}