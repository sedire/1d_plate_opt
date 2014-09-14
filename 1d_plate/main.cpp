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

int main()
{
	cout << "hello\n";

	//omp_set_num_threads( THREAD_NUM );

	//Solver<HPD<N_PRES, GRAD_SIZE> >* solver = new Solver<HPD<N_PRES, GRAD_SIZE> >();

	N_PRES weightJ = 50000.0l;
	N_PRES weightB = 1.0l / 6.0 / 6.0 / 6.0 * 6.0;

	N_PRES J0start =  0.0;
	N_PRES tauStart = 0.0111611;
	N_PRES tauStartExp = 0.0311281;

	N_PRES J0start_1 =  0.02;
	N_PRES tauStart_1 = 0.00516662;
	N_PRES tauStartExp_1 = 0.00915119;

	N_PRES J0start_2 =  0.02;
	N_PRES tauStart_2 = 0.00281418;
	N_PRES tauStartExp_2 = 0.0126425;

	N_PRES J0start_3 =  0.0;
	N_PRES tauStart_3 = 0.00391129;
	N_PRES tauStartExp_3 = 0.0151232;

	N_PRES ByStart = 0.0;

	//N_PRES dJ = 0.0001;
	//N_PRES dTau = 0.000001;
	////double xx[GRAD_SIZE_FULL] = {0.25037, 0.11209, 1e-005, 0.0554859, 0.00446793, 1e-005, 0.285602, 0.187283, 1e-005, -0.135677, 0.312661, 1e-005};
	//double xx[GRAD_SIZE_FULL] = {0.01, 0.0048, 0.0048, 0.01, 0.0048, 0.0048, 0.01, 0.0048, 0.0048, 0.01, 0.0048, 0.0048};
	//double xx1[GRAD_SIZE_FULL] = {0.01, 0.0048, 0.0048 + dTau, 0.01, 0.0048, 0.0048, 0.01, 0.0048, 0.0048, 0.01, 0.0048, 0.0048};
	//double xx2[GRAD_SIZE_FULL] = {0.01, 0.0048, 0.0048 - dTau, 0.01, 0.0048, 0.0048, 0.01, 0.0048, 0.0048, 0.01, 0.0048, 0.0048};
	//double grad[GRAD_SIZE_FULL];
	//for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	//{
	//	grad[i] = 0.0;
	//}
	//N_PRES val1 = calcValTaus( xx, -1 );
	//N_PRES val2 = calcValGradTaus( grad, xx, -1 );
	//cout << " val1 " << val1 << endl;
	//cout << " val2 " << val2 << endl;
	//cout << " diff is " << val1 - val2 << endl;

	//N_PRES f1 = calcValTaus( xx1, -1 );
	//N_PRES f2 = calcValTaus( xx2, -1 );
	//N_PRES deriv = ( f1 - f2 ) / 2.0 / dTau;
	//cout << " finite diff " << deriv << endl;
	//cout << " hpd " << grad[2] << endl;

	//N_PRES*** resArr = new N_PRES**[CHAR_TIME / DELTA_T + 1];		//warning here!!!!
	//for( int i = 0; i < CHAR_TIME / DELTA_T + 1; ++i )
	//{
	//	resArr[i] = new N_PRES*[NODES_Y];

	//	for (int j = 0; j < NODES_Y; ++j)
	//	{
	//		resArr[i][j] = new N_PRES[EQ_NUM];
	//	}
	//}
	Solver<N_PRES>* solver = new Solver<N_PRES>();
	solver->setResArray( 0 );
	solver->setTask( J0start, tauStart, tauStartExp, J0start_3, tauStart_3, tauStartExp_3, ByStart, 10000000, 0.01 );
	solver->setSwitchTime( SWITCH_TIME );

	while( solver->cur_t <= CHAR_TIME )
	{
		cout << solver->cur_t << endl;

		solver->do_step();
		solver->increaseTime();

		solver->dump_check_sol( -1 );
		solver->dump_whole_sol( 4 );
	}

	AdjSolver adjSolver;
	adjSolver.loadParamsFromStruct( solver->saveParamsToStruct() );

	cout << ".........\n";
	cout << "... done!\n";

//////////////////////////////////
	//Solver</*HPD<*/N_PRES/*, GRAD_SIZE>*/ >* solver2 = new Solver</*HPD<*/N_PRES/*, GRAD_SIZE>*/ >();
	//solver2->setTask( J0start, tauStart, tauStartExp, J0start, tauStart, tauStartExp, ByStart, 10000000, 0.01 );
	//time_t tBegin = time( 0 );

	//while( solver2->cur_t <= 0.05 )
	//{
	//	cout << solver2->cur_t << endl;
	//	solver2->do_step();
	//	solver2->cur_t += solver2->dt;
	//	++( solver2->curTimeStep );
	//	solver2->dump_check_sol( -1 );
	//}
	//time_t tEnd = time( 0 );
	//cout << " \n computations are done in " << tEnd - tBegin << endl;
	//cout << ".........\n";
	//cout << "... done!\n";

	cout << " total time: " << solver->totalTime << endl;
	cout << " rgk time: " << solver->rgkTime << endl;
	cout << " matr time: " << solver->matrTime << endl;
	cout << " buildSoln time: " << solver->buildSolnTime << endl;
	cout << "  ortho time: " << solver->orthoTime << endl;
	cout << " ortho time from orthoBuilder: " << solver->getOrthoBTime() << endl;

///////////////////////////////////////

	//Matrix<N_PRES, GRAD_SIZE_FULL, 1> params;
	//params << J0start, tauStart, tauStartExp, 
	//	J0start_1, tauStart_1, tauStartExp_1,
	//	J0start_2, tauStart_2, tauStartExp_2,
	//	J0start_3, tauStart_3, tauStartExp_3;

	//Optimizer<HPD<N_PRES, GRAD_SIZE> > optimizer( solver, weightJ, weightB, CHAR_TIME );
	//optimizeASA_Taus<HPD<N_PRES, GRAD_SIZE> >( params );

	delete solver;
	/*cout << " deleting the solution array now...\n";
	for( int i = 0; i < CHAR_TIME / DELTA_T + 1; ++i )
	{
		for (int j = 0; j < NODES_Y; ++j)
		{
			delete[] resArr[i][j];
		}
		delete[] resArr[i];
	}
	delete[] resArr;*/

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}