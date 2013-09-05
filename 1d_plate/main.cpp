#include <iostream>
#include <stdio.h>
#include "Solver.h"
#include "Plate.h"
#include "Optimizer.h"
#include "hyperDual.h"
#include <omp.h>

using std::cout;
using std::endl;

int main()
{
	cout << "hello\n";

	omp_set_num_threads( 2 );

	//Plate<complex<N_PRES> >* plate = new Plate<complex<N_PRES> >();
	Plate<HPD<N_PRES, GRAD_SIZE> >* plate = new Plate<HPD<N_PRES, GRAD_SIZE> >();

	plate->loadVals( 102970000000.0, 7550000000.0, 0.3, 1594.0, 0.0021, 39000.0, 0.1524 );

	Solver<HPD<N_PRES, GRAD_SIZE> >* solver = new Solver<HPD<N_PRES, GRAD_SIZE> >();
	solver->loadPlate( plate );

	N_PRES weight = 1.0l / 6.0 / 6.0 / 6.0;

	N_PRES J0start =  44000.0 / J0_SCALE;
	N_PRES tauStart = 0.0048;
	N_PRES ByStart = 0.1 / BY0_SCALE;

	Optimizer<HPD<N_PRES, GRAD_SIZE> > optimizer( solver, weight, 0.05 );
	optimizer.optimize( J0start, tauStart, ByStart );

	free( solver );
	free( plate );

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}