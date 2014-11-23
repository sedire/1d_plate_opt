#include <iostream>
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

N_PRES* GlobalResArrays = new N_PRES[1];//[SCEN_NUMBER * ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM];
N_PRES* GlobalResDtArrays = new N_PRES[1];//[SCEN_NUMBER * ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM];
N_PRES* GlobalResAdjArrays = new N_PRES[1];//[SCEN_NUMBER * ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM];

int main()
{
	cout << "hello\n";

	omp_set_num_threads( THREAD_NUM );

	time_t adjTime = 0;
	time_t hpdTime = 0;

	//N_PRES J0start =  0.0215698;
	//N_PRES tauStart = 0.0111486;
	//N_PRES tauStartExp = 0.0198796;

	//N_PRES J0start_1 =  1.0;
	//N_PRES tauStart_1 = 0.00515053;
	//N_PRES tauStartExp_1 = 0.00259337;

	//N_PRES J0start_2 =  0.00991937;
	//N_PRES tauStart_2 = 0.00682416;
	//N_PRES tauStartExp_2 = 0.00001;

	//N_PRES J0start_3 =  1.0;
	//N_PRES tauStart_3 = 0.00370372;
	//N_PRES tauStartExp_3 = 0.00295792;

	N_PRES J0start =  0.0215701;//0.01;
	N_PRES tauStart = 0.0111491;
	N_PRES tauStartExp = 0.0198698;

	N_PRES J0start_1 =  1.0;
	N_PRES tauStart_1 = 0.0051509;
	N_PRES tauStartExp_1 = 0.00259362;

	N_PRES J0start_2 =  0.0100736;
	N_PRES tauStart_2 = 0.00745042;
	N_PRES tauStartExp_2 = 0.00001;

	N_PRES J0start_3 =  1.0;
	N_PRES tauStart_3 = 0.00370382;
	N_PRES tauStartExp_3 = 0.00295845;

	N_PRES ByStart = 1.0;

//////////////////////////////////
	Solver</*HPD<*/N_PRES/*, GRAD_SIZE>*/ >* solver = new Solver</*HPD<*/N_PRES/*, GRAD_SIZE>*/ >();
	solver->setTask( J0start, tauStart, tauStartExp, J0start_2, tauStart_2, tauStartExp_2, ByStart, stress_whole, 100.0, GlobalTauP2 );
	time_t tBegin = time( 0 );

	while( solver->cur_t <= CHAR_TIME )
	{
		cout << solver->cur_t << endl;
		solver->do_step();
		//solver->dump_check_sol( -1 );
		solver->dump_whole_sol( 4 );

		solver->increaseTime();
	}
	time_t tEnd = time( 0 );
	cout << " \n computations are done in " << tEnd - tBegin << endl;
	cout << ".........\n";
	cout << "... done!\n";


	cout << " adjTime " << adjTime << endl;
	cout << " hpdTIme " << hpdTime << endl;

///////////////////////////////////////

	Matrix<N_PRES, GRAD_SIZE_FULL, 1> params;
	params << J0start, tauStart, tauStartExp, 
		J0start_1, tauStart_1, tauStartExp_1,
		J0start_2, tauStart_2, tauStartExp_2,
		J0start_3, tauStart_3, tauStartExp_3;

	//double x[GRAD_SIZE_FULL] = { 0.0267403, 0.0110724, 0.0213401, 
	//							1.0, 0.00514927, 0.00271984,
	//							0.00991937, 0.00682416, 0.00001,
	//							1.0, 0.00371071, 0.00313839 };
	//double x[GRAD_SIZE_FULL] = { 0.01, 0.0048, 0.0048, 
	//							0.01, 0.0048, 0.0048, 
	//							0.01, 0.0048, 0.0048, 
	//							0.01, 0.0048, 0.0048 };

	//double gAdj[GRAD_SIZE_FULL];
	//double g[GRAD_SIZE_FULL];

	//double valAdj = calcValGradTausAdjSolid( gAdj, x, 0 );
	//double val = calcValGradTaus( g, x, 0 );

	//cout << " ----------\n";
	//cout << " == " << val << " " << valAdj << endl;
	//for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	//{
	//	cout << " :: " << g[i] << " " << gAdj[i] << " " << fabs( ( g[i] - gAdj[i] ) / g[i] ) * 100.0 << " % " << endl;
	//}
	//cout << " ----------\n";

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