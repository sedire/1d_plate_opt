#include <iostream>
#include <stdio.h>
#include "Solver.h"
#include "Plate.h"
#include "Optimizer.h"
#include "hyperDual.h"
#include "hyperDual2.h"
#include <omp.h>

using std::cout;
using std::endl;

int main()
{
	cout << "hello\n";

	HPD2<N_PRES, 4> xx1;
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
	return 0;


	omp_set_num_threads( 2 );

	//Plate<complex<N_PRES> >* plate = new Plate<complex<N_PRES> >();
	Plate<HPD<N_PRES, GRAD_SIZE> >* plate = new Plate<HPD<N_PRES, GRAD_SIZE> >();

	plate->loadVals( 102970000000.0, 7550000000.0, 0.3, 1594.0, 0.0021, 39000.0, 0.1524 );

	Solver<HPD<N_PRES, GRAD_SIZE> >* solver = new Solver<HPD<N_PRES, GRAD_SIZE> >();
	solver->loadPlate( plate );

	N_PRES weight = 1.0l / 6.0 / 6.0 / 6.0;

	N_PRES J0start =  0.44553;
	N_PRES tauStart = 44.9132;
	N_PRES ByStart = 11.0971;

	Optimizer<HPD<N_PRES, GRAD_SIZE> > optimizer( solver, weight, 0.05 );
	optimizer.optimize( J0start, tauStart, ByStart );

	delete solver;
	delete plate;

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}