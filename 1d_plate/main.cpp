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
using std::ofstream;

N_PRES* GlobalResArrays = new N_PRES[1];//[SCEN_NUMBER * ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM];
N_PRES* GlobalResDtArrays = new N_PRES[1];//[SCEN_NUMBER * ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM];
N_PRES* GlobalResAdjArrays = new N_PRES[1];//[SCEN_NUMBER * ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM];

int main()
{
	cout << "hello\n";

	omp_set_num_threads( THREAD_NUM );

	time_t adjTime = 0;
	time_t hpdTime = 0;

	N_PRES J0start =  0.01;
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
	N_PRES tauStartExp_3 = 0.0048;

	//N_PRES J0start =  0.0215701;
	//N_PRES tauStart = 0.0111491;
	//N_PRES tauStartExp = 0.0198698;

	//N_PRES J0start_1 =  1.0;
	//N_PRES tauStart_1 = 0.0051509;
	//N_PRES tauStartExp_1 = 0.00259362;

	//N_PRES J0start_2 =  0.0100736;
	//N_PRES tauStart_2 = 0.00745042;
	//N_PRES tauStartExp_2 = 0.00001;

	//N_PRES J0start_3 =  1.0;
	//N_PRES tauStart_3 = 0.00370382;
	//N_PRES tauStartExp_3 = 0.00295845;

	N_PRES ByStart = 1.0;

//////////////////////////////////
//	Solver</*HPD<*/N_PRES/*, GRAD_SIZE>*/ >* solver = new Solver</*HPD<*/N_PRES/*, GRAD_SIZE>*/ >();
//	solver->setTask( J0start, tauStart, tauStartExp, J0start_2, tauStart_2, tauStartExp_2, ByStart, stress_whole, 100.0, GlobalTauP2 );
//	time_t tBegin = time( 0 );
//
//	while( solver->cur_t <= CHAR_TIME )
//	{
//		cout << solver->cur_t << endl;
//		solver->do_step();
//		//solver->dump_check_sol( -1 );
//		solver->dump_whole_sol( 4 );
//
//		solver->increaseTime();
//	}
//	time_t tEnd = time( 0 );
//	cout << " \n computations are done in " << tEnd - tBegin << endl;
//	cout << ".........\n";
//	cout << "... done!\n";
//
//
//	cout << " adjTime " << adjTime << endl;
//	cout << " hpdTIme " << hpdTime << endl;
//
/////////////////////////////////////////

	Matrix<N_PRES, GRAD_SIZE_FULL, 1> params;
	params << J0start, tauStart, tauStartExp, 
		J0start_1, tauStart_1, tauStartExp_1,
		J0start_2, tauStart_2, tauStartExp_2,
		J0start_3, tauStart_3, tauStartExp_3;

	/*double x[GRAD_SIZE_FULL] = { 0.0267403, 0.0110724, 0.0213401, 
								1.0, 0.00514927, 0.00271984,
								0.00991937, 0.00682416, 0.00001,
								1.0, 0.00371071, 0.00313839 };
	double x[GRAD_SIZE_FULL] = { 0.01, 0.0048, 0.0048, 
								0.01, 0.0048, 0.0048, 
								0.01, 0.0048, 0.0048, 
								0.01, 0.0048, 0.0048 };

	double gAdj[GRAD_SIZE_FULL];
	double g[GRAD_SIZE_FULL];

	double valAdj = calcValGradTausAdjSolid( gAdj, x, 0 );
	double val = calcValGradTaus( g, x, 0 );

	cout << " ----------\n";
	cout << " == " << val << " " << valAdj << endl;
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		cout << " :: " << g[i] << " " << gAdj[i] << " " << fabs( ( g[i] - gAdj[i] ) / g[i] ) * 100.0 << " % " << endl;
	}
	cout << " ----------\n";*/

	//CHECKING DERIVATIVE COMPUTATION:
	/*HPD<N_PRES, GRAD_SIZE_SECOND> J0 = J0start;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauSin = tauStart;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauExp = tauStartExp;
	HPD<N_PRES, GRAD_SIZE_SECOND> J0_2 = J0start_2;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauSin_2 = tauStart_2;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauExp_2 = tauStartExp_2;

	J0.elems[1] = 1.0l;
	tauSin.elems[2] = 1.0l;
	tauExp.elems[3] = 1.0l;
	J0_2.elems[4] = 1.0l;
	tauSin_2.elems[5] = 1.0l;
	tauExp_2.elems[6] = 1.0l;

	Solver<HPD<N_PRES, GRAD_SIZE_SECOND> > solverHPD;
	solverHPD.setTask( J0, tauSin, tauExp, J0_2, tauSin_2, tauExp_2, ByStart, stress_centered, GlobalP02, GlobalTauP2 );

	HPD<N_PRES, GRAD_SIZE_SECOND> valHPD;
	HPD<N_PRES, GRAD_SIZE_SECOND> funcValHPD = 0.0l;
	while( solverHPD.cur_t <= CHAR_TIME )
	{
		cout << " t= " << solverHPD.cur_t << endl;
		valHPD = solverHPD.do_step();
		funcValHPD += valHPD * valHPD;
		solverHPD.increaseTime();
	}

	ofstream of1( "HPDderiv.txt" );
	for( int i = 0; i < GRAD_SIZE_SECOND + 1; ++i )
	{
		of1 << funcValHPD.elems[i] << endl;
	}
	of1.close();*/

	/*Solver<N_PRES> solverReg;
	N_PRES dJ = 1e-6;
	N_PRES dTau = 1e-8;
	
	N_PRES val = 0.0;
	N_PRES funcVal1 = 0.0;
	N_PRES funcVal2 = 0.0;

	solverReg.setTask( J0start, tauStart, tauStartExp, J0start_2, tauStart_2, tauStartExp_2 + dTau, ByStart, stress_centered, GlobalP02, GlobalTauP2 );
	while( solverReg.cur_t <= CHAR_TIME )
	{
		cout << " t= " << solverReg.cur_t << endl;
		val = solverReg.do_step();
		funcVal1 += val * val;
		solverReg.increaseTime();
	}

	solverReg.setTask( J0start, tauStart, tauStartExp, J0start_2, tauStart_2, tauStartExp_2 - dTau, ByStart, stress_centered, GlobalP02, GlobalTauP2 );
	while( solverReg.cur_t <= CHAR_TIME )
	{
		cout << " t= " << solverReg.cur_t << endl;
		val = solverReg.do_step();
		funcVal2 += val * val;
		solverReg.increaseTime();
	}

	cout << " fin diff deriv " << ( funcVal1 - funcVal2 ) / ( 2.0 * dTau ) << endl;*/

	optimizeASA_Taus<HPD<N_PRES, GRAD_SIZE> >( params );

	cout << "\n -- Deleting the solution arrays now...\n";

	delete[] GlobalResArrays;
	delete[] GlobalResDtArrays;
	delete[] GlobalResAdjArrays;

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}