#include <iostream>
#include <stdio.h>
#include "Solver.h"
#include "Plate.h"
#include "Optimizer.h"
//#include "hyperDual.h"
//#include "hyperDual2.h"
#include <omp.h>

using std::cout;
using std::endl;

int main()
{
	cout << "hello\n";

	/*HPD2<N_PRES, 4> xx1;
	HPD2<N_PRES, 4> xx2;
	HPD2<N_PRES, 4> xx3;
	HPD2<N_PRES, 4> xx4;

	HPD2<N_PRES, 4> yy;
	
	xx1.elems[0] = 5.01;
	xx1.elems[1] = 1.0;
	xx2.elems[0] = 23.22;
	xx2.elems[2] = 1.0;
	xx3.elems[0] = 4.36;
	xx3.elems[3] = 1.0;
	xx4.elems[0] = 1.6;
	xx4.elems[4] = 1.0;

	yy = xx1 * xx1 + xx2 * xx2 + xx3 * xx3 * xx3 + 2.0l * ( xx1 + xx2 ) + xx4 * ( xx3 - xx2 ) * ( xx3 - xx2 ) * ( xx3 - xx2 ) * ( xx2 - xx1 ) + 20.0l * xx4 * ( xx1 + xx3 ) * ( xx1 + xx3 ) + xx4 * xx4 * xx4 * xx4
		+ 19.0l / xx3 / xx3 / xx3 / xx3 / xx3 - ( xx1 + xx2 + xx3 ) * ( xx2 - 1.0l ) / ( xx4 + xx2 ) / ( xx4 + xx2 )
		+ xx1 * xx1 * ( xx2 - xx3 ) * sin( xx1 + xx2 * cos( exp( xx3 / ( xx4 + xx2 ) / sqrt( xx1 ) ) ) - sqrt( xx2 + xx3 ) );

	cout << yy << endl;

	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;*/

	omp_set_num_threads( 2 );

	//plate->loadVals( 102970000000.0, 7550000000.0, 0.3, 1594.0, 0.0021, 39000.0, 0.1524 );

	//Solver<HPD<N_PRES, GRAD_SIZE> >* solver = new Solver<HPD<N_PRES, GRAD_SIZE> >();

	N_PRES weightJ = 50000.0l;
	N_PRES weightB = 1.0l / 6.0 / 6.0 / 6.0 * 6.0;

	N_PRES J0start =  0.01;
	N_PRES tauStart = 0.0048;
	N_PRES ByStart = 1;

//////////////////////////////////
	Solver<N_PRES>* solver2 = new Solver<N_PRES>();
	solver2->setTask( J0start, tauStart, ByStart );
	time_t tBegin = time( 0 );
	while( solver2->cur_t <= 0.05 )
	{
		cout << solver2->cur_t << endl;
		for( int i = 0; i < 1; ++i )
		{
			solver2->do_step();
			solver2->cur_t += solver2->dt;
			++( solver2->curTimeStep );
		}
		solver2->dump_check_sol( -1 );
	}
	time_t tEnd = time( 0 );
	cout << " \n computations are done in " << tEnd - tBegin << endl;
	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
///////////////////////////////////////

	Matrix<N_PRES, GRAD_SIZE, 1> params;
	params << J0start, tauStart, ByStart;

	//Optimizer<HPD<N_PRES, GRAD_SIZE> > optimizer( solver, weightJ, weightB, CHAR_TIME );
	//optimizer.optimizeASA( params );

	//delete solver;

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}