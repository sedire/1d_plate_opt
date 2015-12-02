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
#include "ParamSack.h"

using std::cout;
using std::endl;
using std::ifstream;

N_PRES* GlobalResArrays = new N_PRES[SCEN_NUMBER * ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM];
N_PRES* GlobalResDtArrays = new N_PRES[SCEN_NUMBER * ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM];
N_PRES* GlobalResAdjArrays = new N_PRES[SCEN_NUMBER * ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM];

ParamSack GlobalParams;

int main()
{
	cout << "hello\n";
	if( GlobalParams.loadFromFile( "params.txt" ) == 0 )
	{
		omp_set_num_threads( THREAD_NUM );

		////////////////////////////////////////
		//Solver<N_PRES>* solver = new Solver<N_PRES>();
		//solver->setTask( GlobalParams.getCurrentType(), GlobalParams.getCurrentParams( 0 ), GlobalParams.getBy0(), 
		//				GlobalParams.getStressType(), GlobalParams.getStressParams( 0 ) );
		//time_t tBegin = time( 0 );
	
		//while( solver->cur_t <= GlobalParams.getTotalTimeT1() )
		//{
		//	cout << solver->cur_t << endl;
		//	solver->do_step();
		//	solver->dumpSolIntegral();
		//	//solver->dumpCheckSol( -1, 0 );
		//	//solver->dumpWholeSol( 1 );
	
		//	solver->increaseTime();
		//}
		//if( solver->getMaxNewtonIterReached() == 1 )
		//{
		//	cout << "max newton iter reached\n";
		//}
	
		//time_t tEnd = time( 0 );
		//cout << " \n computations are done in " << tEnd - tBegin << endl;
		//cout << ".........\n";
		//cout << "... done!\n";
	//
	//
	//	cout << " adjTime " << adjTime << endl;
	//	cout << " hpdTIme " << hpdTime << endl;

	/////////////////////////////////////////

		/*double alphaStep = 1e-2;
		double beta[GRAD_SIZE_FULL];
		beta[0] = 1e-2;
		beta[1] = 1e-4;
		beta[2] = 1e-4;

		beta[3] = -1e-2;
		beta[4] = 1e-4;
		beta[5] = 1e-4;

		beta[6] = 1e-2;
		beta[7] = 1e-4;
		beta[8] = 1e-4;

		beta[9] = -1e-2;
		beta[10] = 1e-4;
		beta[11] = 1e-4;

		double xStart[GRAD_SIZE_FULL];
		xStart[0] = 0.0250578;
		xStart[1] = 0.0111564;
		xStart[2] = 0.0196501;
		xStart[3] = 1;
		xStart[4] = 0.00515453;
		xStart[5] = 0.00267568;
		xStart[6] = 0.010188;
		xStart[7] = 0.00907185;
		xStart[8] = 1e-005;
		xStart[9] = 1;
		xStart[10] = 0.00371065;
		xStart[11] = 0.0030765;

		for( int i = 0; i < 100; ++i )
		{
			for( int j = 0; j < GRAD_SIZE_FULL; ++j )
			{
				xStart[j] += alphaStep * beta[j];
			}
			double fVal = calcValTaus( xStart, 0 );
			ofstream of1( "beta_test.txt", ofstream::out | ofstream::app );
			of1 << i + 1 << " ; " << fVal << endl;
			of1.close();
		}*/

		//------------->>
		//Matrix<N_PRES, GRAD_SIZE_FULL, 1> params;
		//params << GlobalParams.getCurrentParams1st( 0 ), GlobalParams.getCurrentParams1st( 1 ), GlobalParams.getCurrentParams1st( 2 ), 
		//	GlobalParams.getCurrentParams2nd( 0, 0 ), GlobalParams.getCurrentParams2nd( 0, 1 ), GlobalParams.getCurrentParams2nd( 0, 2 ),
		//	GlobalParams.getCurrentParams2nd( 1, 0 ), GlobalParams.getCurrentParams2nd( 1, 1 ), GlobalParams.getCurrentParams2nd( 1, 2 ),
		//	GlobalParams.getCurrentParams2nd( 2, 0 ), GlobalParams.getCurrentParams2nd( 2, 1 ), GlobalParams.getCurrentParams2nd( 2, 2 );

		//double x[GRAD_SIZE_FULL] = { 0.0267403, 0.0110724, 0.0213401, 
		//							1.0, 0.00514927, 0.00271984,
		//							0.00991937, 0.00682416, 0.00001,
		//							1.0, 0.00371071, 0.00313839 };
		//double x[GRAD_SIZE_FULL] = { 0.01, 0.0048, 0.0048, 
		//							0.01, 0.0048, 0.0048, 
		//							0.01, 0.0048, 0.0048, 
		//							0.01, 0.0048, 0.0048 };
		//double dtt = 0.0000001;
		//double x1[GRAD_SIZE_FULL] = { 0.01, 0.0048 - dtt, 0.0048, 
		//							0.01, 0.0048, 0.0048, 
		//							0.01, 0.0048, 0.0048, 
		//							0.01, 0.0048, 0.0048 };
		//double x2[GRAD_SIZE_FULL] = { 0.01, 0.0048 + dtt, 0.0048, 
		//							0.01, 0.0048, 0.0048, 
		//							0.01, 0.0048, 0.0048, 
		//							0.01, 0.0048, 0.0048 };

		//double gAdj[GRAD_SIZE_FULL];
		//double g[GRAD_SIZE_FULL];

		//double valAdj = calcValGradTausAdj( gAdj, x, 0 );
		//cout << " adj comput done\n";
		//double val = calcValGradTaus( g, x, 0 );

		//ifstream iff( "HPDgrad.txt" );
		//if( iff.is_open() )
		//{
		//	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
		//	{
		//		iff >> g[i];
		//	}
		//	iff.close();
		//}

		//cout << " ----------\n";
		//cout << " == " << val /*<< " " << valAdj*/ << endl;
		//for( int i = 0; i < GRAD_SIZE_FULL; ++i )
		//{
		//	cout << " :: " << g[i] /*<< " " << gAdj[i] << " " << fabs( ( g[i] - gAdj[i] ) / g[i] ) * 100.0 << " % "*/ << endl;
		//}
		//cout << " ----------\n";

		//double val1 = calcValTaus( x1, 0 );
		//double val2 = calcValTaus( x2, 0 );
		//cout << " finite diff " << ( val2 - val1 ) / 2.0 / dtt << endl;

		optimizeASA_Taus();
		//optimizeASAPiece();

		/*int resArrSize = ( ( int )( CHAR_TIME / DELTA_T  + 1 ) ) * NODES_Y * EQ_NUM;
		double gAdj[26];
		double x[26] = { 2.0, -2.0, 
							-0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
							-0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
							-0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, };
		double dJJ = 1e-7;
		double x1[26] = { 2.0, -2.0, 
					-0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
					-0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
					-0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		double x2[26] = { 2.0, -2.0, 
					-0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
					-0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
					-0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		double x11[26];
		double x22[26];

		for( int i = 5; i < 26; ++i )
		{
			for( int j = 0; j < 26; ++j )
			{
				x11[j] = x1[j];
				x22[j] = x2[j];
			}

			x11[i] -= dJJ;
			x22[i] += dJJ;

			double val1 = calcValTaus( x11, 0 );
			double val2 = calcValTaus( x22, 0 );
			cout << " finite diff " << ( val2 - val1 ) / 2.0 / dJJ << endl;
			ofstream of ( "finDiffDeriv.txt", ofstream::app );
			of << ( val2 - val1 ) / 2.0 / dJJ << endl;
		}*/
		//ofstream of( "adjDeriv.txt" );
		//for( int i = 0; i < 26; ++i )
		//{
		//	of << gAdj[i] << endl;
		//}
		//of.close();
	}
	else
	{
		cout << "ERROR loading parameters from file. Terminating...\n";
	}

	cout << "\n -- Deleting the solution arrays now...\n";

	delete[] GlobalResArrays;
	delete[] GlobalResDtArrays;
	delete[] GlobalResAdjArrays;

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}