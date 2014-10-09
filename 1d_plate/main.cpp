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

	time_t adjTime = 0;
	time_t hpdTime = 0;

	//omp_set_num_threads( THREAD_NUM );

	//N_PRES weightJ = 50000.0l;
	//N_PRES weightB = 1.0l / 6.0 / 6.0 / 6.0 * 6.0;

	adjTime = time( 0 );

	N_PRES J0start =  0.01;
	N_PRES tauStart = 0.048;
	N_PRES tauStartExp = 0.048;

	N_PRES J0start_1 =  0.02;
	N_PRES tauStart_1 = 0.00516662;
	N_PRES tauStartExp_1 = 0.00915119;

	N_PRES J0start_2 =  0.02;
	N_PRES tauStart_2 = 0.00281418;
	N_PRES tauStartExp_2 = 0.0126425;

	N_PRES J0start_3 =  0.0;
	N_PRES tauStart_3 = 0.00391129;
	N_PRES tauStartExp_3 = 0.0151232;

	N_PRES ByStart = 1.0;

	N_PRES p0 = 10000000.0;
	N_PRES tauP = 0.01;

	N_PRES* resArr = new N_PRES[( CHAR_TIME / DELTA_T + 1 ) * NODES_Y * EQ_NUM];		//warning here!!!!
	for( int i = 0; i < ( CHAR_TIME / DELTA_T + 1 ) * NODES_Y * EQ_NUM; ++i )
	{
		resArr[i] = 0.0;
	}
	N_PRES* resArrDt = new N_PRES[( CHAR_TIME / DELTA_T + 1 ) * NODES_Y * EQ_NUM];		//warning here!!!!
	for( int i = 0; i < ( CHAR_TIME / DELTA_T + 1 ) * NODES_Y * EQ_NUM; ++i )
	{
		resArrDt[i] = 0.0;
	}
	N_PRES* resArrAdj = new N_PRES[( CHAR_TIME / DELTA_T + 1 ) * NODES_Y * EQ_NUM];		//warning here!!!!
	for( int i = 0; i < ( CHAR_TIME / DELTA_T + 1 ) * NODES_Y * EQ_NUM; ++i )
	{
		resArrAdj[i] = 0.0;
	}
	Solver<N_PRES>* solver = new Solver<N_PRES>();
	solver->setResArray( resArr );
	solver->setResArrayDt( resArrDt );
	solver->setTask( J0start, tauStart, tauStartExp, J0start_1, tauStart_1, tauStartExp_1, ByStart, p0, tauP );
	solver->setSwitchTime( SWITCH_TIME );

	while( solver->cur_t <= CHAR_TIME )
	{
		cout << solver->cur_t << endl;

		solver->do_step();

		//solver->dump_check_sol( -1 );
		solver->dumpSolAll( -1 );
		//solver->dump_whole_sol( 4 );

		solver->increaseTime();
	}

	cout << "creating the AdjSolver\n";
	AdjSolver adjSolver;
	adjSolver.loadParamsFromStruct( solver->saveParamsToStruct() );
	adjSolver.setPrimalSolnData( resArr );
	adjSolver.setPrimalDtData( resArrDt );
	adjSolver.setAdjointSolnData( resArrAdj );

	while( adjSolver.getCurTimeStep() >= 0 )
	{
		cout << " time is " << adjSolver.getCurTime() << " " << adjSolver.getCurTimeStep() << endl;
		adjSolver.doStep();
		adjSolver.dumpSol( -1 );
		adjSolver.decreaseTime();

		cout << " new time is " << adjSolver.getCurTime() << " " << adjSolver.getCurTimeStep() << endl;
		//std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	}

	adjTime = time( 0 ) - adjTime;

	N_PRES dy = 0.1524 / ( NODES_Y - 1 );
	N_PRES dt = DELTA_T;
	HPD<N_PRES, GRAD_SIZE_SECOND> sum = 0.0l;

	hpdTime = time( 0 );

	HPD<N_PRES, GRAD_SIZE_SECOND> J0startHPD =  J0start;
	J0startHPD.elems[1] = 1.0;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauStartHPD = tauStart;
	tauStartHPD.elems[2] = 1.0;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauStartExpHPD = tauStartExp;
	tauStartExpHPD.elems[3] = 1.0;

	HPD<N_PRES, GRAD_SIZE_SECOND> J0startHPD_1 =  J0start_1;
	J0startHPD_1.elems[4] = 1.0;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauStartHPD_1 = tauStart_1;
	tauStartHPD_1.elems[5] = 1.0;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauStartExpHPD_1 = tauStartExp_1;
	tauStartExpHPD_1.elems[6] = 1.0;

	HPD<N_PRES, GRAD_SIZE_SECOND> ByStartHPD = ByStart;

	Solver<HPD<N_PRES, GRAD_SIZE_SECOND> >* solverHPD = new Solver<HPD<N_PRES, GRAD_SIZE_SECOND> >();
	solverHPD->setTask( J0startHPD, tauStartHPD, tauStartExpHPD, J0startHPD_1, tauStartHPD_1, tauStartExpHPD_1, ByStartHPD, p0, tauP );
	solverHPD->setSwitchTime( SWITCH_TIME );
	/*while( solverHPD->cur_t <= CHAR_TIME )
	{
		cout << solverHPD->cur_t << endl;
		
		if( solverHPD->cur_t <= CHAR_TIME )
		{
			sum += solverHPD->do_step();
		}
		solverHPD->increaseTime();
	}
	sum *= dt * dy;*/

	hpdTime = time( 0 ) - hpdTime;

	N_PRES dJ = 0.00001;
	N_PRES dTau = 0.0000005;
	N_PRES obj1 = 0.0;
	N_PRES obj2 = 0.0;

	/*solver->setTask( J0start, tauStart, tauStartExp + dTau, J0start_3, tauStart_3, tauStartExp_3, ByStart, p0, tauP );
	while( solver->cur_t <= CHAR_TIME )
	{
		cout << solver->cur_t << endl;

		solver->do_step();
		solver->increaseTime();
	}
	for( int t = 0; t < CHAR_TIME / DELTA_T + 1 - 1; ++t )
	{
		for( int y = 0; y < NODES_Y - 1; ++y )
		{
			obj1 += resArr[ NODES_Y * EQ_NUM * t + y * EQ_NUM + 1] * resArr[ NODES_Y * EQ_NUM * t + y * EQ_NUM + 1] * dy * DELTA_T;
		}
	}

	solver->setTask( J0start, tauStart, tauStartExp - dTau, J0start_3, tauStart_3, tauStartExp_3, ByStart, p0, tauP );
	while( solver->cur_t <= CHAR_TIME )
	{
		cout << solver->cur_t << endl;

		solver->do_step();
		solver->increaseTime();
	}
	for( int t = 0; t < CHAR_TIME / DELTA_T + 1 - 1; ++t )
	{
		for( int y = 0; y < NODES_Y - 1; ++y )
		{
			obj2 += resArr[ NODES_Y * EQ_NUM * t + y * EQ_NUM + 1] * resArr[ NODES_Y * EQ_NUM * t + y * EQ_NUM + 1] * dy * DELTA_T;
		}
	}*/
	
	//cout << " dF / dTauExp by fin diff  " << ( obj1 - obj2 ) / ( 2.0 * dTau ) << endl;
	cout << " dF / dJ0 " <<  adjSolver.calcJ0Deriv() << " " << sum.elems[1] << endl;
	cout << " dF / dJ1 " <<  adjSolver.calcJ1Deriv() << " " << sum.elems[4] << endl;
	cout << " dF / dtauSin0 " <<  adjSolver.calcTauSin0Deriv() << " " << sum.elems[2] << endl;
	cout << " dF / dtauSin1 " <<  adjSolver.calcTauSin1Deriv() << " " << sum.elems[5] << endl;
	cout << " dF / dtauExp0 " <<  adjSolver.calcTauExp0Deriv() << " " << sum.elems[3] << endl;
	cout << " dF / dtauExp1 " <<  adjSolver.calcTauExp1Deriv() << " " << sum.elems[6] << endl;

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
	cout << "ortho time from orthoBuilder: " << solver->getOrthoBTime() << endl;

	cout << " adjTime " << adjTime << endl;
	cout << " hpdTIme " << hpdTime << endl;

///////////////////////////////////////

	//Matrix<N_PRES, GRAD_SIZE_FULL, 1> params;
	//params << J0start, tauStart, tauStartExp, 
	//	J0start_1, tauStart_1, tauStartExp_1,
	//	J0start_2, tauStart_2, tauStartExp_2,
	//	J0start_3, tauStart_3, tauStartExp_3;

	//Optimizer<HPD<N_PRES, GRAD_SIZE> > optimizer( solver, weightJ, weightB, CHAR_TIME );
	//optimizeASA_Taus<HPD<N_PRES, GRAD_SIZE> >( params );

	delete solver;
	cout << "\n -- Deleting the solution arrays now...\n";
	delete[] resArr;
	delete[] resArrDt;
	delete[] resArrAdj;

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}