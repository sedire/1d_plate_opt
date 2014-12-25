#include "HagerOptFuncs.h"

extern ParamSack GlobalParams;

double calcValGradTausAdj( double* g, double* x, long n )
{
	time_t begin = time( 0 );

	double ret = 0.0l;
	N_PRES dt = DELTA_T;
	N_PRES dy = 0.1524 / ( NODES_Y - 1 );
	int resArrSize = ( ( int )( CHAR_TIME / DELTA_T  + 1 ) ) * NODES_Y * EQ_NUM;

	cout << "try to calc at\n";
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		cout << x[i] << endl;
	}
	cout << " ====\n";

	N_PRES J0begin1;
	N_PRES tauBeginSin1;
	N_PRES tauBeginExp1;

	N_PRES J0begin2[SCEN_NUMBER];
	N_PRES tauBeginSin2[SCEN_NUMBER];
	N_PRES tauBeginExp2[SCEN_NUMBER];

	N_PRES B0begin2;

	J0begin1 = x[0];
	tauBeginSin1 = x[1];
	tauBeginExp1 = x[2];
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		J0begin2[scen] = x[( scen + 1 ) * 3];
		tauBeginSin2[scen] = x[( scen + 1 ) * 3 + 1];
		tauBeginExp2[scen] = x[( scen + 1 ) * 3 + 2];
	}
	B0begin2 = 1.0l;

	Solver<N_PRES> solver[SCEN_NUMBER];
	AdjSolver adjSolver[SCEN_NUMBER];

	cout << "\tcalculating func val\n";

	N_PRES funcVal1[SCEN_NUMBER];
	N_PRES funcVal2[SCEN_NUMBER];
	N_PRES mechLoad[SCEN_NUMBER] = { GlobalP01, GlobalP02, GlobalP03 };
	N_PRES mechTaus[SCEN_NUMBER] = { GlobalTauP1, GlobalTauP2, GlobalTauP3 };

#pragma omp parallel for
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		cout << omp_get_thread_num() << endl;
		
		solver[scen].setTask( J0begin1, tauBeginSin1, tauBeginExp1, J0begin2[scen], tauBeginSin2[scen], tauBeginExp2[scen], 
								B0begin2, stressCentered, GlobalParams.getStressParams( scen, 0 ), mechLoad[scen], mechTaus[scen] );
		solver[scen].setResArray( GlobalResArrays + scen * resArrSize );
		solver[scen].setResDtArray( GlobalResDtArrays + scen * resArrSize );

		N_PRES sum = 0.0;
		while( solver[scen].cur_t <= SWITCH_TIME )
		{
			sum += solver[scen].do_step();
			solver[scen].increaseTime();
		}
		funcVal1[scen] = sum * dt * dy / SWITCH_TIME;

		sum = 0.0;
		while( solver[scen].cur_t <= CHAR_TIME )
		{
			sum += solver[scen].do_step();
			solver[scen].increaseTime();
		}
		funcVal2[scen] = sum * dt * dy / ( CHAR_TIME - SWITCH_TIME );

		if( solver[scen].getMaxNewtonIterReached() == 1 )
		{
			cout << " WARNING: scenario " << scen << " solver used too many newton iterations\n";
		}

		adjSolver[scen].loadParamsFromStruct( solver[scen].saveParamsToStruct() );
		adjSolver[scen].setPrimalSolnData( GlobalResArrays + scen * resArrSize );
		adjSolver[scen].setPrimalDtData( GlobalResDtArrays + scen * resArrSize );
		adjSolver[scen].setAdjointSolnData( GlobalResAdjArrays + scen * resArrSize );

		while( adjSolver[scen].getCurTime() >= SWITCH_TIME )
		{
			//cout << " --- " << adjSolver[scen].getCurTime() << endl;
			adjSolver[scen].doStep();
			adjSolver[scen].decreaseTime();
		}
		adjSolver[scen].scalePrevTimeStepSoln();
		while( adjSolver[scen].getCurTime() >= 0.0 )
		{
			adjSolver[scen].doStep();
			adjSolver[scen].decreaseTime();
		}
	}

	if( g != 0 )
	{
		g[0] = ( adjSolver[0].calcJ0DerivS() + adjSolver[1].calcJ0DerivS() + adjSolver[2].calcJ0DerivS() ) / 3.0;

		g[1] = ( adjSolver[0].calcTauSin0DerivS() + adjSolver[1].calcTauSin0DerivS() + adjSolver[2].calcTauSin0DerivS() ) / 3.0;

		g[2] = ( adjSolver[0].calcTauExp0DerivS() + adjSolver[1].calcTauExp0DerivS() + adjSolver[2].calcTauExp0DerivS() ) / 3.0;

		g[3] = adjSolver[0].calcJ1DerivS() / 3.0;
		g[4] = adjSolver[0].calcTauSin1DerivS() / 3.0;
		g[5] = adjSolver[0].calcTauExp1DerivS() / 3.0;

		g[6] = adjSolver[1].calcJ1DerivS() / 3.0;
		g[7] = adjSolver[1].calcTauSin1DerivS() / 3.0;
		g[8] = adjSolver[1].calcTauExp1DerivS() / 3.0;

		g[9] = adjSolver[2].calcJ1DerivS() / 3.0;
		g[10] = adjSolver[2].calcTauSin1DerivS() / 3.0;
		g[11] = adjSolver[2].calcTauExp1DerivS() / 3.0;
	}

	ret = ( funcVal1[0] + funcVal1[1] + funcVal1[2]
			+ funcVal2[0] + funcVal2[1] + funcVal2[2] ) / 3.0l;

	time_t endtime = time( 0 );
	cout << "\tdone in " << endtime - begin << endl;

	return ret;
}

double calcValGradTausAdjSolid( double* g, double* x, long n )	//"Solid" means that there is only one integral in the objective, both stages are considered inside that integral
{
	time_t begin = time( 0 );

	double ret = 0.0l;
	N_PRES dt = DELTA_T;
	N_PRES dy = 0.1524 / ( NODES_Y - 1 );
	int resArrSize = ( ( int )( CHAR_TIME / DELTA_T  + 1 ) ) * NODES_Y * EQ_NUM;

	cout << "try to calc at\n";
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		cout << x[i] << endl;
	}
	cout << " ====\n";

	N_PRES J0begin1;
	N_PRES tauBeginSin1;
	N_PRES tauBeginExp1;

	N_PRES J0begin2[SCEN_NUMBER];
	N_PRES tauBeginSin2[SCEN_NUMBER];
	N_PRES tauBeginExp2[SCEN_NUMBER];

	N_PRES B0begin2;

	J0begin1 = x[0];
	tauBeginSin1 = x[1];
	tauBeginExp1 = x[2];
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		J0begin2[scen] = x[( scen + 1 ) * 3];
		tauBeginSin2[scen] = x[( scen + 1 ) * 3 + 1];
		tauBeginExp2[scen] = x[( scen + 1 ) * 3 + 2];
	}
	B0begin2 = 1.0l;

	Solver<N_PRES> solver[SCEN_NUMBER];
	AdjSolver adjSolver[SCEN_NUMBER];

	cout << "\tcalculating func val\n";

	N_PRES funcVal1[SCEN_NUMBER];
	N_PRES funcVal2[SCEN_NUMBER];
	N_PRES mechLoad[SCEN_NUMBER] = { GlobalP01, GlobalP02, GlobalP03 };
	N_PRES mechTaus[SCEN_NUMBER] = { GlobalTauP1, GlobalTauP2, GlobalTauP3 };

#pragma omp parallel for
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		cout << omp_get_thread_num() << endl;
		
		solver[scen].setTask( J0begin1, tauBeginSin1, tauBeginExp1, J0begin2[scen], tauBeginSin2[scen], tauBeginExp2[scen],
								B0begin2, stressCentered, GlobalParams.getStressParams( scen, 0 ), mechLoad[scen], mechTaus[scen] );
		solver[scen].setResArray( GlobalResArrays + scen * resArrSize );
		solver[scen].setResDtArray( GlobalResDtArrays + scen * resArrSize );

		N_PRES sum = 0.0;
		while( solver[scen].cur_t <= SWITCH_TIME )
		{
			sum += solver[scen].do_step();
			solver[scen].increaseTime();
		}
		funcVal1[scen] = sum * dt * dy;

		sum = 0.0;
		while( solver[scen].cur_t <= CHAR_TIME )
		{
			sum += solver[scen].do_step();
			solver[scen].increaseTime();
		}
		funcVal2[scen] = sum * dt * dy;

		if( solver[scen].getMaxNewtonIterReached() == 1 )
		{
			cout << " WARNING: scenario " << scen << " solver used too many newton iterations\n";
		}

		adjSolver[scen].loadParamsFromStruct( solver[scen].saveParamsToStruct() );
		adjSolver[scen].setPrimalSolnData( GlobalResArrays + scen * resArrSize );
		adjSolver[scen].setPrimalDtData( GlobalResDtArrays + scen * resArrSize );
		adjSolver[scen].setAdjointSolnData( GlobalResAdjArrays + scen * resArrSize );

		while( adjSolver[scen].getCurTime() >= 0.0 )
		{
			//cout << " --- " << adjSolver[scen].getCurTime() << endl;
			adjSolver[scen].doStep();
			adjSolver[scen].decreaseTime();
		}
	}

	if( g != 0 )
	{
		g[0] = ( adjSolver[0].calcJ0Deriv() + adjSolver[1].calcJ0Deriv() + adjSolver[2].calcJ0Deriv() ) / 3.0;

		g[1] = ( adjSolver[0].calcTauSin0Deriv() + adjSolver[1].calcTauSin0Deriv() + adjSolver[2].calcTauSin0Deriv() ) / 3.0;

		g[2] = ( adjSolver[0].calcTauExp0Deriv() + adjSolver[1].calcTauExp0Deriv() + adjSolver[2].calcTauExp0Deriv() ) / 3.0;

		g[3] = adjSolver[0].calcJ1Deriv() / 3.0;
		g[4] = adjSolver[0].calcTauSin1Deriv() / 3.0;
		g[5] = adjSolver[0].calcTauExp1Deriv() / 3.0;

		g[6] = adjSolver[1].calcJ1Deriv() / 3.0;
		g[7] = adjSolver[1].calcTauSin1Deriv() / 3.0;
		g[8] = adjSolver[1].calcTauExp1Deriv() / 3.0;

		g[9] = adjSolver[2].calcJ1Deriv() / 3.0;
		g[10] = adjSolver[2].calcTauSin1Deriv() / 3.0;
		g[11] = adjSolver[2].calcTauExp1Deriv() / 3.0;
	}

	ret = ( funcVal1[0] + funcVal1[1] + funcVal1[2]
			+ funcVal2[0] + funcVal2[1] + funcVal2[2] ) / 3.0l;

	time_t endtime = time( 0 );
	cout << "\tdone in " << endtime - begin << endl;

	return ret;
}

double calcValGradTaus( double* g, double* x, long n )
{
	time_t begin = time( 0 );

	double ret = 0.0l;
	N_PRES dt = DELTA_T;
	N_PRES dy = 0.1524 / ( NODES_Y - 1 );

	cout << "try to calc at\n";
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		cout << x[i] << endl;
	}
	cout << " ====\n";

	HPD<N_PRES, GRAD_SIZE_SECOND> J0begin1;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauBeginSin1;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauBeginExp1;

	HPD<N_PRES, GRAD_SIZE_SECOND> J0begin2[SCEN_NUMBER];
	HPD<N_PRES, GRAD_SIZE_SECOND> tauBeginSin2[SCEN_NUMBER];
	HPD<N_PRES, GRAD_SIZE_SECOND> tauBeginExp2[SCEN_NUMBER];

	HPD<N_PRES, GRAD_SIZE_SECOND> B0begin2;

	J0begin1.elems[0] = x[0];
	J0begin1.elems[1] = 1.0l;
	J0begin1.elems[2] = 0.0l;
	J0begin1.elems[3] = 0.0l;
	J0begin1.elems[4] = 0.0l;
	J0begin1.elems[5] = 0.0l;
	J0begin1.elems[6] = 0.0l;

	tauBeginSin1.elems[0] = x[1];
	tauBeginSin1.elems[1] = 0.0l;
	tauBeginSin1.elems[2] = 1.0l;
	tauBeginSin1.elems[3] = 0.0l;
	tauBeginSin1.elems[4] = 0.0l;
	tauBeginSin1.elems[5] = 0.0l;
	tauBeginSin1.elems[6] = 0.0l;

	tauBeginExp1.elems[0] = x[2];
	tauBeginExp1.elems[1] = 0.0l;
	tauBeginExp1.elems[2] = 0.0l;
	tauBeginExp1.elems[3] = 1.0l;
	tauBeginExp1.elems[4] = 0.0l;
	tauBeginExp1.elems[5] = 0.0l;
	tauBeginExp1.elems[6] = 0.0l;

	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		J0begin2[scen].elems[0] = x[( scen + 1 ) * 3];
		J0begin2[scen].elems[1] = 0.0l;
		J0begin2[scen].elems[2] = 0.0l;
		J0begin2[scen].elems[3] = 0.0l;
		J0begin2[scen].elems[4] = 1.0l;
		J0begin2[scen].elems[5] = 0.0l;
		J0begin2[scen].elems[6] = 0.0l;

		tauBeginSin2[scen].elems[0] = x[( scen + 1 ) * 3 + 1];
		tauBeginSin2[scen].elems[1] = 0.0l;
		tauBeginSin2[scen].elems[2] = 0.0l;
		tauBeginSin2[scen].elems[3] = 0.0l;
		tauBeginSin2[scen].elems[4] = 0.0l;
		tauBeginSin2[scen].elems[5] = 1.0l;
		tauBeginSin2[scen].elems[6] = 0.0l;

		tauBeginExp2[scen].elems[0] = x[( scen + 1 ) * 3 + 2];
		tauBeginExp2[scen].elems[1] = 0.0l;
		tauBeginExp2[scen].elems[2] = 0.0l;
		tauBeginExp2[scen].elems[3] = 0.0l;
		tauBeginExp2[scen].elems[4] = 0.0l;
		tauBeginExp2[scen].elems[5] = 0.0l;
		tauBeginExp2[scen].elems[6] = 1.0l;
	}

	B0begin2 = 1.0l;

	Solver<HPD<N_PRES, GRAD_SIZE_SECOND> > solver_second[SCEN_NUMBER];

	cout << "\tcalculating func and grad val\n";

	double charTime = CHAR_TIME;

	HPD<N_PRES, GRAD_SIZE_SECOND> funcVal1[SCEN_NUMBER];
	HPD<N_PRES, GRAD_SIZE_SECOND> funcVal2[SCEN_NUMBER];
	N_PRES mechLoad[SCEN_NUMBER] = { GlobalP01, GlobalP02, GlobalP03 };
	N_PRES mechTaus[SCEN_NUMBER] = { GlobalTauP1, GlobalTauP2, GlobalTauP3 };

#pragma omp parallel for
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		solver_second[scen].setTask( J0begin1, tauBeginSin1, tauBeginExp1, J0begin2[scen], tauBeginSin2[scen], tauBeginExp2[scen],
										B0begin2, stressCentered, GlobalParams.getStressParams( scen, 0 ), mechLoad[scen], mechTaus[scen] );

		HPD<N_PRES, GRAD_SIZE_SECOND> sum = 0.0;
		funcVal1[scen] = 0.0l;
		while( solver_second[scen].cur_t <= SWITCH_TIME )
		{
			//sum += solver_second[scen].do_step();
			sum = solver_second[scen].do_step();
			funcVal1[scen] += sum * sum;

			solver_second[scen].increaseTime();
		}
		//funcVal1[scen] = sum * dt * dy / SWITCH_TIME;
		funcVal1[scen] /= SWITCH_TIME;

		funcVal2[scen] = 0.0l;
		sum = 0.0;
		while( solver_second[scen].cur_t <= charTime )
		{
			//sum += solver_second[scen].do_step();
			sum = solver_second[scen].do_step();
			funcVal2[scen] += sum * sum; 

			solver_second[scen].increaseTime();
		}
		//funcVal2[scen] = sum * dt * dy / ( charTime - SWITCH_TIME );
		funcVal2[scen] /= ( charTime - SWITCH_TIME );

		if( solver_second[scen].getMaxNewtonIterReached() == 1 )
		{
			cout << " WARNING: scenario " << scen << " solver used too many newton iterations\n";
		}
	}

	N_PRES Weight = J_WEIGHT;
	if( g != 0 )
	{
		g[0] = ( funcVal1[0].elems[1] + funcVal1[1].elems[1] + funcVal1[2].elems[1]
				+ funcVal2[0].elems[1] + funcVal2[1].elems[1] + funcVal2[2].elems[1] ) / 3.0;

		g[1] = ( funcVal1[0].elems[2] + funcVal1[1].elems[2] + funcVal1[2].elems[2]
				+ funcVal2[0].elems[2] + funcVal2[1].elems[2] + funcVal2[2].elems[2] ) / 3.0;

		g[2] = ( funcVal1[0].elems[3] + funcVal1[1].elems[3] + funcVal1[2].elems[3]
				+ funcVal2[0].elems[3] + funcVal2[1].elems[3] + funcVal2[2].elems[3] ) / 3.0;

		g[3] = funcVal2[0].elems[4] / 3.0;
		g[4] = funcVal2[0].elems[5] / 3.0;
		g[5] = funcVal2[0].elems[6] / 3.0;

		g[6] = funcVal2[1].elems[4] / 3.0;
		g[7] = funcVal2[1].elems[5] / 3.0;
		g[8] = funcVal2[1].elems[6] / 3.0;

		g[9] = funcVal2[2].elems[4] / 3.0;
		g[10] = funcVal2[2].elems[5] / 3.0;
		g[11] = funcVal2[2].elems[6] / 3.0;
	}
	ret = ( funcVal1[0].real() + funcVal1[1].real() + funcVal1[2].real()
			+ funcVal2[0].real() + funcVal2[1].real() + funcVal2[2].real() ) / 3.0l;

	cout << " the grad is:\n";
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		cout << "\t" << g[i] << endl;
	}
	cout << "-----------";

	time_t endtime = time( 0 );
	cout << "\tdone in " << endtime - begin << endl;

	return ret;
}

double calcValGradTausDet( double* g, double* x, long n )
{
	time_t begin = time( 0 );

	double ret = 0.0l;
	N_PRES dt = DELTA_T;
	N_PRES dy = 0.1524 / ( NODES_Y - 1 );

	cout << "try to calc at\n";
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		cout << x[i] << endl;
	}
	cout << " ====\n";

	HPD<N_PRES, GRAD_SIZE_SECOND> J0begin1;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauBeginSin1;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauBeginExp1;

	HPD<N_PRES, GRAD_SIZE_SECOND> J0begin2;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauBeginSin2;
	HPD<N_PRES, GRAD_SIZE_SECOND> tauBeginExp2;

	HPD<N_PRES, GRAD_SIZE_SECOND> B0begin2;

	J0begin1.elems[0] = x[0];
	J0begin1.elems[1] = 1.0l;
	J0begin1.elems[2] = 0.0l;
	J0begin1.elems[3] = 0.0l;
	J0begin1.elems[4] = 0.0l;
	J0begin1.elems[5] = 0.0l;
	J0begin1.elems[6] = 0.0l;

	tauBeginSin1.elems[0] = x[1];
	tauBeginSin1.elems[1] = 0.0l;
	tauBeginSin1.elems[2] = 1.0l;
	tauBeginSin1.elems[3] = 0.0l;
	tauBeginSin1.elems[4] = 0.0l;
	tauBeginSin1.elems[5] = 0.0l;
	tauBeginSin1.elems[6] = 0.0l;

	tauBeginExp1.elems[0] = x[2];
	tauBeginExp1.elems[1] = 0.0l;
	tauBeginExp1.elems[2] = 0.0l;
	tauBeginExp1.elems[3] = 1.0l;
	tauBeginExp1.elems[4] = 0.0l;
	tauBeginExp1.elems[5] = 0.0l;
	tauBeginExp1.elems[6] = 0.0l;

	J0begin2.elems[0] = x[( 1 + 1 ) * 3];
	J0begin2.elems[1] = 0.0l;
	J0begin2.elems[2] = 0.0l;
	J0begin2.elems[3] = 0.0l;
	J0begin2.elems[4] = 1.0l;
	J0begin2.elems[5] = 0.0l;
	J0begin2.elems[6] = 0.0l;

	tauBeginSin2.elems[0] = x[( 1 + 1 ) * 3 + 1];
	tauBeginSin2.elems[1] = 0.0l;
	tauBeginSin2.elems[2] = 0.0l;
	tauBeginSin2.elems[3] = 0.0l;
	tauBeginSin2.elems[4] = 0.0l;
	tauBeginSin2.elems[5] = 1.0l;
	tauBeginSin2.elems[6] = 0.0l;

	tauBeginExp2.elems[0] = x[( 1 + 1 ) * 3 + 2];
	tauBeginExp2.elems[1] = 0.0l;
	tauBeginExp2.elems[2] = 0.0l;
	tauBeginExp2.elems[3] = 0.0l;
	tauBeginExp2.elems[4] = 0.0l;
	tauBeginExp2.elems[5] = 0.0l;
	tauBeginExp2.elems[6] = 1.0l;

	B0begin2 = 1.0l;

	Solver<HPD<N_PRES, GRAD_SIZE_SECOND> > solver_second;

	cout << "\tcalculating func and grad val\n";

	double charTime = CHAR_TIME;

	HPD<N_PRES, GRAD_SIZE_SECOND> funcVal1;
	HPD<N_PRES, GRAD_SIZE_SECOND> funcVal2;
	N_PRES mechRad = GlobalParams.getStressParams( 1, 0 );
	N_PRES mechLoad = GlobalP02;
	N_PRES mechTaus = GlobalTauP2;

	solver_second.setTask( J0begin1, tauBeginSin1, tauBeginExp1, J0begin2, tauBeginSin2, tauBeginExp2,
									B0begin2, stressCentered, mechRad, mechLoad, mechTaus );

	HPD<N_PRES, GRAD_SIZE_SECOND> sum = 0.0;
	funcVal1 = 0.0l;
	while( solver_second.cur_t <= SWITCH_TIME )
	{
		sum = solver_second.do_step();
		funcVal1 += sum * sum;

		solver_second.increaseTime();
	}
	funcVal1 /= SWITCH_TIME;

	funcVal2 = 0.0l;
	sum = 0.0;
	while( solver_second.cur_t <= charTime )
	{
		sum = solver_second.do_step();
		funcVal2 += sum * sum; 

		solver_second.increaseTime();
	}
	funcVal2 /= ( charTime - SWITCH_TIME );

	if( solver_second.getMaxNewtonIterReached() == 1 )
	{
		cout << " WARNING: solver used too many newton iterations\n";
	}

	N_PRES Weight = J_WEIGHT;
	if( g != 0 )
	{
		g[0] = funcVal1.elems[1] + funcVal2.elems[1];

		g[1] = funcVal1.elems[2] + funcVal2.elems[2];

		g[2] = funcVal1.elems[3] + funcVal2.elems[3];

		g[3] = 0.0l;
		g[4] = 0.0l;
		g[5] = 0.0l;

		g[6] = funcVal2.elems[4];
		g[7] = funcVal2.elems[5];
		g[8] = funcVal2.elems[6];

		g[9] = 0.0l;
		g[10] = 0.0l;
		g[11] = 0.0l;
	}
	ret = funcVal1.real() + funcVal2.real();

	cout << " the grad is:\n";
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		cout << "\t" << g[i] << endl;
	}
	cout << "-----------";

	time_t endtime = time( 0 );
	cout << "\tdone in " << endtime - begin << endl;

	return ret;
}

double calcValTaus( double* x, long n )
{
	double ret = 0.0l;
	time_t begin = time( 0 );
	N_PRES dt = DELTA_T;
	N_PRES dy = 0.1524 / ( NODES_Y - 1 );

	cout << "try to calc at\n";
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		cout << x[i] << endl;
	}
	cout << " ====\n";

	ofstream of( "ineq_check.txt", ofstream::app );
	of << "calcValTaus\n";
	double lb = 0.001;
	if( x[1] < lb || x[4] < lb || x[7] < lb || x[10] < lb )
	{
		of << "ineq tauSin >= " << lb << " violated\n\n";
		for( int i = 0; i < GRAD_SIZE_FULL; ++i )
		{
			of << x[i] << endl;
		}
		of << "----------------------\n\n";
	}
	of.close();

	N_PRES J0begin;
	N_PRES tauBeginSin;
	N_PRES tauBeginExp;

	N_PRES J0begin2[SCEN_NUMBER];
	N_PRES tauBeginSin2[SCEN_NUMBER];
	N_PRES tauBeginExp2[SCEN_NUMBER];

	N_PRES B0begin;

	J0begin = x[0];
	tauBeginSin = x[1];
	tauBeginExp = x[2];

	for( int i = 0; i < SCEN_NUMBER; ++i )
	{
		J0begin2[i] = x[( i + 1 ) * 3];
		tauBeginSin2[i] = x[( i + 1 ) * 3 + 1];
		tauBeginExp2[i] = x[( i + 1 ) * 3 + 2];
	}
	B0begin = 1.0l;

	Solver<N_PRES> solver[SCEN_NUMBER];

	cout << "\tcalculating func val\n";

	double charTime = CHAR_TIME;

	N_PRES funcVal1[SCEN_NUMBER];
	N_PRES funcVal2[SCEN_NUMBER];
	N_PRES mechLoad[SCEN_NUMBER] = { GlobalP01, GlobalP02, GlobalP03 };
	N_PRES mechTaus[SCEN_NUMBER] = { GlobalTauP1, GlobalTauP2, GlobalTauP3 };

#pragma omp parallel for
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		funcVal1[scen] = 0.0l;
		funcVal2[scen] = 0.0l;
		solver[scen].setTask( J0begin, tauBeginSin, tauBeginExp, J0begin2[scen], tauBeginSin2[scen], tauBeginExp2[scen],
								B0begin, stressCentered, GlobalParams.getStressParams( scen, 0 ), mechLoad[scen], mechTaus[scen] );
		N_PRES sum = 0.0;

		while( solver[scen].cur_t <= SWITCH_TIME )
		{
			//sum += solver[scen].do_step();
			sum = solver[scen].do_step();
			funcVal1[scen] += sum * sum; 

			solver[scen].increaseTime(); 

			//solver_second.dump_check_sol( -1 );
		}
		//funcVal1[scen] = sum * dt * dy / SWITCH_TIME;
		funcVal1[scen] /= SWITCH_TIME;

		sum = 0.0;
		while( solver[scen].cur_t <= charTime )
		{
			//sum += solver[scen].do_step();
			sum = solver[scen].do_step();
			funcVal2[scen] += sum * sum;

			solver[scen].increaseTime();

			//solver_second.dump_check_sol( -1 );
		}
		//funcVal2[scen] = sum * dt * dy / ( charTime - SWITCH_TIME );
		funcVal2[scen] /= ( charTime - SWITCH_TIME );

		if( solver[scen].getMaxNewtonIterReached() == 1 )
		{
			cout << " WARNING: scenario " << scen << " solver used too many newton iterations\n";
		}
	}

	cout << "\tfunc val done\n";
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		cout << "\t1st; scen " << scen << " " << funcVal1[scen] << endl;
		cout << "\t2nd; scen " << scen << " " << funcVal2[scen] << endl;
	}
	cout << " -------------\n";

	N_PRES Weight = J_WEIGHT;
	ret = ( funcVal1[0] + funcVal1[1] + funcVal1[2]
		+ funcVal2[0] + funcVal2[1] + funcVal2[2] ) / 3.0l;

	time_t endtime = time( 0 );
	cout << "\tdone in " << endtime - begin << endl;

	cout << funcVal1[0] + funcVal2[0] << endl;
	cout << funcVal1[1] + funcVal2[1] << endl;
	cout << funcVal1[2] + funcVal2[2] << endl;

	return ret;
}

double calcValTausDet( double* x, long n )
{
	double ret = 0.0l;
	time_t begin = time( 0 );
	N_PRES dt = DELTA_T;
	N_PRES dy = 0.1524 / ( NODES_Y - 1 );

	cout << "try to calc at\n";
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		cout << x[i] << endl;
	}
	cout << " ====\n";

	N_PRES J0begin;
	N_PRES tauBeginSin;
	N_PRES tauBeginExp;

	N_PRES J0begin2;
	N_PRES tauBeginSin2;
	N_PRES tauBeginExp2;

	N_PRES B0begin;

	J0begin = x[0];
	tauBeginSin = x[1];
	tauBeginExp = x[2];

	J0begin2 = x[( 1 + 1 ) * 3];
	tauBeginSin2 = x[( 1 + 1 ) * 3 + 1];
	tauBeginExp2 = x[( 1 + 1 ) * 3 + 2];

	B0begin = 1.0l;

	Solver<N_PRES> solver;

	cout << "\tcalculating func val\n";

	double charTime = CHAR_TIME;

	N_PRES funcVal1;
	N_PRES funcVal2;
	N_PRES mechRad = GlobalParams.getStressParams( 1, 0 );
	N_PRES mechLoad = GlobalP02;
	N_PRES mechTaus = GlobalTauP2;

	funcVal1 = 0.0l;
	funcVal2 = 0.0l;
	solver.setTask( J0begin, tauBeginSin, tauBeginExp, J0begin2, tauBeginSin2, tauBeginExp2,
							B0begin, stressCentered, mechRad, mechLoad, mechTaus );
	N_PRES sum = 0.0;

	while( solver.cur_t <= SWITCH_TIME )
	{
		sum = solver.do_step();
		funcVal1 += sum * sum; 

		solver.increaseTime(); 
	}
	funcVal1 /= SWITCH_TIME;

	sum = 0.0;
	while( solver.cur_t <= charTime )
	{
		sum = solver.do_step();
		funcVal2 += sum * sum;

		solver.increaseTime();
	}
	funcVal2 /= ( charTime - SWITCH_TIME );

	if( solver.getMaxNewtonIterReached() == 1 )
	{
		cout << " WARNING: solver used too many newton iterations\n";
	}

	N_PRES Weight = J_WEIGHT;
	ret = funcVal1 + funcVal2;

	time_t endtime = time( 0 );
	cout << "\tdone in " << endtime - begin << endl;

	return ret;
}

double calc1stOrdOptInfoASA_Taus( asa_objective* asa )
{
	cout << "\tcalc 1st order CG_DES Both\n";
	double ret = calcValGradTaus( asa->g, asa->x, asa->n );
	return ret;
}

double calcValASA_Taus( asa_objective* asa )
{
	cout << "\tcalc 1st order CG_DES Val\n";
	return calcValTaus( asa->x, asa->n );
}

void calcGradASA_Taus( asa_objective* asa )
{
	cout << "\tcalc 1st order CG_DES Grad\n";
	calcValGradTaus( asa->g, asa->x, asa->n );
}
