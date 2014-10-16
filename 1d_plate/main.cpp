#include <iostream>
#include <stdio.h>
#include "Solver.h"
#include "Plate.h"
#include "Optimizer.h"
//#include "hyperDual.h"
//#include "hyperDual2.h"
#include <omp.h>
#include "HagerOptFuncs.h"

using std::cout;
using std::endl;

int main()
{
	cout << "hello\n";

	omp_set_num_threads( THREAD_NUM );

	Solver<HPD<N_PRES, GRAD_SIZE> >* solver = new Solver<HPD<N_PRES, GRAD_SIZE> >();

	//N_PRES weightJ = 50000.0l;
	//N_PRES weightB = 1.0l / 6.0 / 6.0 / 6.0 * 6.0;

	N_PRES J0start =  0.0215698;
	N_PRES tauStart = 0.0111486;
	N_PRES tauStartExp = 0.0198796;

	N_PRES J0start_1 =  1.0;
	N_PRES tauStart_1 = 0.00515053;
	N_PRES tauStartExp_1 = 0.00259337;

	N_PRES J0start_2 =  0.00991937;
	N_PRES tauStart_2 = 0.00682416;
	N_PRES tauStartExp_2 = 0.00001;

	N_PRES J0start_3 =  1.0;
	N_PRES tauStart_3 = 0.00370372;
	N_PRES tauStartExp_3 = 0.00295792;

	N_PRES ByStart = 1.0;

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

	//cout << " total time: " << solver2->totalTime << endl;
	//cout << " rgk time: " << solver2->rgkTime << endl;
	//cout << " matr time: " << solver2->matrTime << endl;
	//cout << " buildSoln time: " << solver2->buildSolnTime << endl;
	//cout << "  ortho time: " << solver2->orthoTime << endl;
	//cout << " ortho time from orthoBuilder: " << solver2->getOrthoBTime() << endl;

	//std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	//return 0;
///////////////////////////////////////

	Matrix<N_PRES, GRAD_SIZE_FULL, 1> params;
	params << J0start, tauStart, tauStartExp, 
		J0start_1, tauStart_1, tauStartExp_1,
		J0start_2, tauStart_2, tauStartExp_2,
		J0start_3, tauStart_3, tauStartExp_3;

	double x[GRAD_SIZE_FULL] = { 0.0213812, 0.0111541, 0.0198321, 
								1.0, 0.00515865, 0.0025922,
								0.00991937, 0.00682416, 0.00001,
								1.0, 0.00370119, 0.00295973 };
	calcValTaus( x, 0 );

	//Optimizer<HPD<N_PRES, GRAD_SIZE> > optimizer( solver, weightJ, weightB, CHAR_TIME );
	//optimizeASA_Taus<HPD<N_PRES, GRAD_SIZE> >( params );

	delete solver;

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}