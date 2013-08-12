#include <iostream>
#include <stdio.h>
#include "Solver.h"
#include "Plate.h"
#include "Optimizer.h"
#include <omp.h>

using std::cout;
using std::endl;

int main()
{
	cout << "hello\n";

	omp_set_num_threads( 2 );

	cout << sizeof( float ) << " " << sizeof( double ) << " " << sizeof( long double ) << endl;

	Plate<complex<N_PRES> >* plate = new Plate<complex<N_PRES> >();
	plate->loadVals( 102970000000.0, 7550000000.0, 0.3, 1594.0, 0.0021, 39000.0, 0.1524 );

	Solver<complex<N_PRES> >* solver = new Solver<complex<N_PRES> >();
	solver->loadPlate( plate );

	N_PRES weight = 1.0l / 6.0 / 6.0 / 6.0;
	N_PRES J0h = 0.00001;
	N_PRES tauh = 0.00001;
	N_PRES J0start =  44000.0 / J0_SCALE;
	N_PRES tauStart = 0.0048;
	Optimizer<complex<N_PRES> > optimizer( solver, weight, J0h, tauh, 0.05 );
	optimizer.optimize( J0start, tauStart );

	free( solver );
	free( plate );

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}