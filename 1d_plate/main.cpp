#include <iostream>
#include <stdio.h>
#include "Solver.h"
#include "Plate.h"
#include "Optimizer.h"
//#include "hyperDual.h"
#include "HagerOptFuncs.h"
//#include "hyperDual2.h"
#include <omp.h>

using std::cout;
using std::endl;

int main()
{
	cout << "hello\n";

	//omp_set_num_threads( THREAD_NUM );

	Solver<HPD<N_PRES, GRAD_SIZE> >* solver = new Solver<HPD<N_PRES, GRAD_SIZE> >();

	N_PRES weightJ = 50000.0l;
	N_PRES weightB = 1.0l / 6.0 / 6.0 / 6.0 * 6.0;

	N_PRES J0start =  0.01;
	N_PRES tauStart = 0.0048;
	N_PRES tauStartExp = 0.0048;

	N_PRES J0start_1 =  0.01;
	N_PRES tauStart_1 = 0.0048;
	N_PRES tauStartExp_1 = 0.0048;

	N_PRES J0start_2 =  0.01;
	N_PRES tauStart_2 = 0.0048;
	N_PRES tauStartExp_2 = 0.0048;

	N_PRES J0start_3 =  0.01;
	N_PRES tauStart_3 = 0.0048;
	N_PRES tauStartExp_3 = 0.0048;

	N_PRES ByStart = 1;

	HPD<N_PRES, 2> a = 1.0;
	HPD<N_PRES, 2> b = 2.0;
	HPD<N_PRES, 2> c = 3.0;
	HPD<N_PRES, 2> d = 4.0;
	HPD<N_PRES, 2> e = 5.0;

	a.elems[1] = 2.0;
	a.elems[2] = 3.0;
	b.elems[1] = 4.0;
	b.elems[2] = 5.0;
	c.elems[1] = 6.0;
	c.elems[2] = 7.0;
	d.elems[1] = 8.0;
	d.elems[2] = 9.0;
	e.elems[1] = 10.0;
	e.elems[2] = 11.0;

	cout << a << endl;
	cout << b << endl;
	cout << c << endl;
	cout << d << endl;
	cout << e << endl;
	cout << " --------------\n";

	a *= b;
	cout << a << endl;
	a *= 3.0l;
	cout << a << endl;
	a /= d;
	cout << a << endl;
	a /= 5.0l; 
	cout << a << endl;
	a /= a;
	cout << a << endl;
	a = b;
	cout << a << endl;
	a *= a;
	cout << a << endl;

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


	/*Solver<N_PRES> solver1;
	Solver<HPD<N_PRES, GRAD_SIZE_FIRST> > solver2a;
	Solver<HPD<N_PRES, GRAD_SIZE_SECOND> > solver2b;

	solver1.setTask( 0.02, 0.008, 0.005, 0.03, 0.008, 0.004548, 1.0, 10000000, 0.01 );
	solver1.setMechLoad( 7500000.0, 0.008 );

	HPD<N_PRES, GRAD_SIZE_FIRST> J0_1a = 0.02;
	J0_1a.elems[1] = 1.0;
	HPD<N_PRES, GRAD_SIZE_FIRST> tauSin_1a = 0.008;
	tauSin_1a.elems[2] = 1.0;
	HPD<N_PRES, GRAD_SIZE_FIRST> tauExp_1a = 0.005;
	tauExp_1a.elems[3] = 1.0;

	HPD<N_PRES, GRAD_SIZE_FIRST> J0_2a = 0.03;
	HPD<N_PRES, GRAD_SIZE_FIRST> tauSin_2a = 0.008;
	HPD<N_PRES, GRAD_SIZE_FIRST> tauExp_2a = 0.004548;

	solver2a.setTask( J0_1a, tauSin_1a, tauExp_1a, J0_2a, tauSin_2a, tauExp_2a, 1.0, 10000000, 0.01 );
	solver2a.setMechLoad( 7500000.0, 0.008 );

	HPD<N_PRES, GRAD_SIZE_SECOND> J0_1b = 0.02;
	J0_1b.elems[1] = 1.0;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauSin_1b = 0.008;
	tauSin_1b.elems[2] = 1.0;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauExp_1b = 0.005;
	tauExp_1b.elems[3] = 1.0;

	HPD<N_PRES, GRAD_SIZE_SECOND> J0_2b = 0.03;
	J0_2b.elems[4] = 1.0;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauSin_2b = 0.008;
	tauSin_2b.elems[5] = 1.0;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauExp_2b = 0.004548;
	tauExp_2b.elems[6] = 1.0;

	solver2b.setTask( J0_1b, tauSin_1b, tauExp_1b, J0_2b, tauSin_2b, tauExp_2b, 1.0, 10000000, 0.01 );
	solver2b.setMechLoad( 7500000.0, 0.008 );

	N_PRES val1 = 0.0;
	HPD<N_PRES, GRAD_SIZE_FIRST> val2a = 0.0;
	HPD<N_PRES, GRAD_SIZE_SECOND> val2b = 0.0;

	N_PRES funcVal1a = 0.0;
	N_PRES funcVal1b = 0.0;
	HPD<N_PRES, GRAD_SIZE_FIRST> funcVal2a = 0.0;
	HPD<N_PRES, GRAD_SIZE_SECOND> funcVal2b = 0.0;

	while( solver1.cur_t <= SWITCH_TIME )
	{
		cout << " == " << solver1.cur_t << endl;
		val1 = solver1.do_step();
		val2a = solver2a.do_step();

		funcVal1a += val1 * val1;
		funcVal2a += val2a * val2a;

		solver1.increaseTime();
		solver2a.increaseTime();
	}
	funcVal1a /= SWITCH_TIME;
	funcVal2a /= SWITCH_TIME;

	while( solver1.cur_t <= CHAR_TIME )
	{
		cout << " == " << solver1.cur_t << endl;
		val1 = solver1.do_step();

		funcVal1b += val1 * val1;

		solver1.increaseTime();
	}
	funcVal1b /= ( CHAR_TIME - SWITCH_TIME );

	while( solver2b.cur_t <= SWITCH_TIME )
	{
		cout << " == " << solver2b.cur_t << endl;
		solver2b.do_step();

		solver2b.increaseTime();
	}
	while( solver2b.cur_t <= CHAR_TIME )
	{
		cout << " == " << solver2b.cur_t << endl;
		val2b = solver2b.do_step();
		funcVal2b += val2b * val2b;

		solver2b.increaseTime();
	}
	funcVal2b /= ( CHAR_TIME - SWITCH_TIME );

	cout << " val1 " << funcVal1a + funcVal1b << endl;
	cout << " val2 " << funcVal2a.real() + funcVal2b.real() << endl;
	cout << " diff is " << funcVal1a + funcVal1b - funcVal2a.real() - funcVal2b.real() << endl;*/


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

	//Optimizer<HPD<N_PRES, GRAD_SIZE> > optimizer( solver, weightJ, weightB, CHAR_TIME );
	//optimizeASA_Taus<HPD<N_PRES, GRAD_SIZE> >( params );

	delete solver;

	cout << ".........\n";
	cout << "... done!\n";
	std::cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
	return 0;
}