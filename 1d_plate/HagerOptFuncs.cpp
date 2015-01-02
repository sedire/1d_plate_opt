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

	vector<vector<N_PRES> > currentParams;
	currentParams.resize( GlobalParams.getNumberOfScenarios(), vector<N_PRES>( GlobalParams.getNumberOfCurrentParams(), 0.0 ) );
	int firstStageParamNum = -1;
	if( GlobalParams.getCurrentType() == currentExpSin )
	{
		firstStageParamNum = 3;
	}
	else if( GlobalParams.getCurrentType() == currentPieceLin )
	{
		N_PRES dT = CHAR_TIME / GlobalParams.getNumberOfCurrentParams();
		firstStageParamNum = (int)floor( SWITCH_TIME / dT );
	}
	else
	{
		cout << "ERROR: we do not support other current types yet\n";
	}

	for( int scen = 0; scen < GlobalParams.getNumberOfScenarios(); ++scen )
	{
		for( int i = 0; i < firstStageParamNum; ++i )
		{
			currentParams[scen][i] = x[i];
		}
		for( int i = firstStageParamNum; i < GlobalParams.getNumberOfCurrentParams(); ++i )
		{
			currentParams[scen][i] = x[scen * ( GlobalParams.getNumberOfCurrentParams() - firstStageParamNum ) + i];
		}
	}

	N_PRES B0begin2 = GlobalParams.getBy0();

	Solver<N_PRES> solver[SCEN_NUMBER];
	AdjSolver adjSolver[SCEN_NUMBER];

	cout << "\tcalculating func val\n";

	N_PRES funcVal1[SCEN_NUMBER];
	N_PRES funcVal2[SCEN_NUMBER];

#pragma omp parallel for
	for( int scen = 0; scen < GlobalParams.getNumberOfScenarios(); ++scen )
	{
		cout << omp_get_thread_num() << endl;
		
		solver[scen].setTask( GlobalParams.getCurrentType(), currentParams[scen], 
								B0begin2, GlobalParams.getStressType(), GlobalParams.getStressParams( scen ) );
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
		if( GlobalParams.getCurrentType() == currentExpSin )
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
		else if( GlobalParams.getCurrentType() == currentPieceLin )
		{
			N_PRES dT = CHAR_TIME / GlobalParams.getNumberOfCurrentParams();
			int firstStageParamNum = (int)floor( SWITCH_TIME / dT );
			int secondStageParamNum = GlobalParams.getNumberOfCurrentParams() - firstStageParamNum;
			for( int i = 0; i < firstStageParamNum; ++i )
			{
				N_PRES sum = 0.0;
				for( int scen = 0; scen < GlobalParams.getNumberOfScenarios(); ++scen )
				{
					sum += adjSolver[scen].calcPieceLinDeriv( i );
				}
				g[i] = sum / GlobalParams.getNumberOfScenarios();
			}
			for( int scen = 0; scen < GlobalParams.getNumberOfScenarios(); ++scen )
			{
				for( int i = firstStageParamNum; i < GlobalParams.getNumberOfCurrentParams(); ++i )
				{
					g[i + scen * secondStageParamNum] = adjSolver[scen].calcPieceLinDeriv( i ) / GlobalParams.getNumberOfScenarios();
				}
			}
		}
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

	vector<vector<N_PRES> > currentParams;
	currentParams.resize( GlobalParams.getNumberOfScenarios(), vector<N_PRES>( GlobalParams.getNumberOfCurrentParams(), 0.0 ) );
	int firstStageParamNum = -1;
	if( GlobalParams.getCurrentType() == currentExpSin )
	{
		firstStageParamNum = 3;
	}
	else if( GlobalParams.getCurrentType() == currentPieceLin )
	{
		N_PRES dT = CHAR_TIME / GlobalParams.getNumberOfCurrentParams();
		firstStageParamNum = (int)( SWITCH_TIME / dT );
	}
	else
	{
		cout << "ERROR: we do not support other current types yet\n";
	}

	for( int scen = 0; scen < GlobalParams.getNumberOfScenarios(); ++scen )
	{
		for( int i = 0; i < firstStageParamNum; ++i )
		{
			currentParams[scen][i] = x[i];
		}
		for( int i = firstStageParamNum; i < GlobalParams.getNumberOfCurrentParams(); ++i )
		{
			currentParams[scen][i] = x[scen * ( GlobalParams.getNumberOfCurrentParams() - firstStageParamNum ) + i];
		}
	}

	N_PRES B0begin2 = GlobalParams.getBy0();

	Solver<N_PRES> solver[SCEN_NUMBER];
	AdjSolver adjSolver[SCEN_NUMBER];

	cout << "\tcalculating func val\n";

	N_PRES funcVal1[SCEN_NUMBER];
	N_PRES funcVal2[SCEN_NUMBER];

#pragma omp parallel for
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		cout << omp_get_thread_num() << endl;
		
		solver[scen].setTask( GlobalParams.getCurrentType(), currentParams[scen], 
								B0begin2, GlobalParams.getStressType(), GlobalParams.getStressParams( scen ) );
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

	vector<vector<HPD<N_PRES, GRAD_SIZE_SECOND> > > currentParams;
	currentParams.resize( GlobalParams.getNumberOfScenarios(), vector<HPD<N_PRES, GRAD_SIZE_SECOND> >( GlobalParams.getNumberOfCurrentParams(), 0.0l ) );
	int firstStageParamNum = -1;
	if( GlobalParams.getCurrentType() == currentExpSin )
	{
		firstStageParamNum = 3;
	}
	else if( GlobalParams.getCurrentType() == currentPieceLin )
	{
		N_PRES dT = CHAR_TIME / GlobalParams.getNumberOfCurrentParams();
		firstStageParamNum = (int)( SWITCH_TIME / dT );
	}
	else
	{
		cout << "ERROR: we do not support other current types yet\n";
	}

	for( int scen = 0; scen < GlobalParams.getNumberOfScenarios(); ++scen )
	{
		for( int i = 0; i < firstStageParamNum; ++i )
		{
			currentParams[scen][i] = x[i];
		}
		for( int i = firstStageParamNum; i < GlobalParams.getNumberOfCurrentParams(); ++i )
		{
			currentParams[scen][i] = x[scen * ( GlobalParams.getNumberOfCurrentParams() - firstStageParamNum ) + i];
		}
	}

	if( GlobalParams.getCurrentType() == currentExpSin )
	{
		for( int scen = 0; scen < GlobalParams.getNumberOfScenarios(); ++scen )
		{
			for( int i = 0; i < GlobalParams.getNumberOfCurrentParams(); ++i )
			{
				currentParams[scen][i].elems[i + 1] = 1.0l; 
			}
		}
	}
	else
	{
		cout << "ERROR: this type of current is not supported yet for HPD method\n";
	}

	HPD<N_PRES, GRAD_SIZE_SECOND> B0begin2 = GlobalParams.getBy0();

	Solver<HPD<N_PRES, GRAD_SIZE_SECOND> > solver[SCEN_NUMBER];

	cout << "\tcalculating func and grad val\n";

	double charTime = CHAR_TIME;

	HPD<N_PRES, GRAD_SIZE_SECOND> funcVal1[SCEN_NUMBER];
	HPD<N_PRES, GRAD_SIZE_SECOND> funcVal2[SCEN_NUMBER];

#pragma omp parallel for
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		solver[scen].setTask( GlobalParams.getCurrentType(), currentParams[scen],
										B0begin2, GlobalParams.getStressType(), GlobalParams.getStressParams( scen ) );

		HPD<N_PRES, GRAD_SIZE_SECOND> sum = 0.0l;
		funcVal1[scen] = 0.0l;
		while( solver[scen].cur_t <= SWITCH_TIME )
		{
			//sum += solver[scen].do_step();
			sum = solver[scen].do_step();
			funcVal1[scen] += sum * sum;

			solver[scen].increaseTime();
		}
		//funcVal1[scen] = sum * dt * dy / SWITCH_TIME;
		funcVal1[scen] /= SWITCH_TIME;

		funcVal2[scen] = 0.0l;
		sum = 0.0;
		while( solver[scen].cur_t <= charTime )
		{
			//sum += solver[scen].do_step();
			sum = solver[scen].do_step();
			funcVal2[scen] += sum * sum; 

			solver[scen].increaseTime();
		}
		//funcVal2[scen] = sum * dt * dy / ( charTime - SWITCH_TIME );
		funcVal2[scen] /= ( charTime - SWITCH_TIME );

		if( solver[scen].getMaxNewtonIterReached() == 1 )
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

	if( GlobalParams.getNumberOfScenarios() != 1 )
	{
		cout << "ERROR: the number of scenarios is != 1 in the deterministic case\n";
		return -1.0;
	}

	vector<vector<HPD<N_PRES, GRAD_SIZE_SECOND> > > currentParams;
	currentParams.resize( GlobalParams.getNumberOfScenarios(), vector<HPD<N_PRES, GRAD_SIZE_SECOND> >( GlobalParams.getNumberOfCurrentParams(), 0.0l ) );
	int firstStageParamNum = -1;
	if( GlobalParams.getCurrentType() == currentExpSin )
	{
		firstStageParamNum = 3;
	}
	else if( GlobalParams.getCurrentType() == currentPieceLin )
	{
		N_PRES dT = CHAR_TIME / GlobalParams.getNumberOfCurrentParams();
		firstStageParamNum = (int)( SWITCH_TIME / dT );
	}
	else
	{
		cout << "ERROR: we do not support other current types yet\n";
	}

	for( int scen = 0; scen < GlobalParams.getNumberOfScenarios(); ++scen )
	{
		for( int i = 0; i < firstStageParamNum; ++i )
		{
			currentParams[scen][i] = x[i];
		}
		for( int i = firstStageParamNum; i < GlobalParams.getNumberOfCurrentParams(); ++i )
		{
			currentParams[scen][i] = x[scen * ( GlobalParams.getNumberOfCurrentParams() - firstStageParamNum ) + i];
		}
	}

	if( GlobalParams.getCurrentType() == currentExpSin )
	{
		for( int scen = 0; scen < GlobalParams.getNumberOfScenarios(); ++scen )
		{
			for( int i = 0; i < GlobalParams.getNumberOfCurrentParams(); ++i )
			{
				currentParams[scen][i].elems[i + 1] = 1.0l; 
			}
		}
	}
	else
	{
		cout << "ERROR: this type of current is not supported yet for HPD method\n";
	}

	HPD<N_PRES, GRAD_SIZE_SECOND> B0begin2 = GlobalParams.getBy0();

	Solver<HPD<N_PRES, GRAD_SIZE_SECOND> > solver;

	cout << "\tcalculating func and grad val\n";

	double charTime = CHAR_TIME;

	HPD<N_PRES, GRAD_SIZE_SECOND> funcVal1;
	HPD<N_PRES, GRAD_SIZE_SECOND> funcVal2;

	solver.setTask( GlobalParams.getCurrentType(), currentParams[0],
									B0begin2, GlobalParams.getStressType(), GlobalParams.getStressParams( 0 ) );

	HPD<N_PRES, GRAD_SIZE_SECOND> sum = 0.0;
	funcVal1 = 0.0l;
	while( solver.cur_t <= SWITCH_TIME )
	{
		sum = solver.do_step();
		funcVal1 += sum * sum;

		solver.increaseTime();
	}
	funcVal1 /= SWITCH_TIME;

	funcVal2 = 0.0l;
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
	
	vector<vector<N_PRES> > currentParams;
	currentParams.resize( GlobalParams.getNumberOfScenarios(), vector<N_PRES>( GlobalParams.getNumberOfCurrentParams(), 0.0 ) );
	int firstStageParamNum = -1;
	if( GlobalParams.getCurrentType() == currentExpSin )
	{
		firstStageParamNum = 3;
	}
	else if( GlobalParams.getCurrentType() == currentPieceLin )
	{
		N_PRES dT = CHAR_TIME / GlobalParams.getNumberOfCurrentParams();
		firstStageParamNum = (int)( SWITCH_TIME / dT );
	}
	else
	{
		cout << "ERROR: we do not support other current types yet\n";
	}

	for( int scen = 0; scen < GlobalParams.getNumberOfScenarios(); ++scen )
	{
		for( int i = 0; i < firstStageParamNum; ++i )
		{
			currentParams[scen][i] = x[i];
		}
		for( int i = firstStageParamNum; i < GlobalParams.getNumberOfCurrentParams(); ++i )
		{
			currentParams[scen][i] = x[scen * ( GlobalParams.getNumberOfCurrentParams() - firstStageParamNum ) + i];
		}
	}

	N_PRES B0begin = GlobalParams.getBy0();

	Solver<N_PRES> solver[SCEN_NUMBER];

	cout << "\tcalculating func val\n";

	double charTime = CHAR_TIME;

	N_PRES funcVal1[SCEN_NUMBER];
	N_PRES funcVal2[SCEN_NUMBER];

#pragma omp parallel for
	for( int scen = 0; scen < SCEN_NUMBER; ++scen )
	{
		funcVal1[scen] = 0.0l;
		funcVal2[scen] = 0.0l;
		solver[scen].setTask( GlobalParams.getCurrentType(), currentParams[scen],
								B0begin, GlobalParams.getStressType(), GlobalParams.getStressParams( scen ) );
		N_PRES sum = 0.0;

		while( solver[scen].cur_t <= SWITCH_TIME )
		{
			sum += solver[scen].do_step();
			//sum = solver[scen].do_step();
			//funcVal1[scen] += sum * sum; 

			solver[scen].increaseTime(); 

			//solver_second.dump_check_sol( -1 );
		}
		funcVal1[scen] = sum * dt * dy / SWITCH_TIME;
		//funcVal1[scen] /= SWITCH_TIME;

		sum = 0.0;
		while( solver[scen].cur_t <= charTime )
		{
			sum += solver[scen].do_step();
			//sum = solver[scen].do_step();
			//funcVal2[scen] += sum * sum;

			solver[scen].increaseTime();

			//solver_second.dump_check_sol( -1 );
		}
		funcVal2[scen] = sum * dt * dy / ( charTime - SWITCH_TIME );
		//funcVal2[scen] /= ( charTime - SWITCH_TIME );

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

	if( GlobalParams.getNumberOfScenarios() != 1 )
	{
		cout << "ERROR: the number of scenarios is != 1 in the deterministic case\n";
		return -1.0;
	}

	vector<vector<N_PRES> > currentParams;
	currentParams.resize( GlobalParams.getNumberOfScenarios(), vector<N_PRES>( GlobalParams.getNumberOfCurrentParams(), 0.0 ) );
	int firstStageParamNum = -1;
	if( GlobalParams.getCurrentType() == currentExpSin )
	{
		firstStageParamNum = 3;
	}
	else if( GlobalParams.getCurrentType() == currentPieceLin )
	{
		N_PRES dT = CHAR_TIME / GlobalParams.getNumberOfCurrentParams();
		firstStageParamNum = (int)( SWITCH_TIME / dT );
	}
	else
	{
		cout << "ERROR: we do not support other current types yet\n";
	}

	for( int scen = 0; scen < GlobalParams.getNumberOfScenarios(); ++scen )
	{
		for( int i = 0; i < firstStageParamNum; ++i )
		{
			currentParams[scen][i] = x[i];
		}
		for( int i = firstStageParamNum; i < GlobalParams.getNumberOfCurrentParams(); ++i )
		{
			currentParams[scen][i] = x[scen * ( GlobalParams.getNumberOfCurrentParams() - firstStageParamNum ) + i];
		}
	}

	N_PRES B0begin = GlobalParams.getBy0();

	Solver<N_PRES> solver;

	cout << "\tcalculating func val\n";

	double charTime = CHAR_TIME;

	N_PRES funcVal1;
	N_PRES funcVal2;

	funcVal1 = 0.0l;
	funcVal2 = 0.0l;
	solver.setTask( GlobalParams.getCurrentType(), currentParams[0],
							B0begin, GlobalParams.getStressType(), GlobalParams.getStressParams( 0 ) );
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

double calc1stOrdOptInfoASAPiece( asa_objective* asa )
{
	cout << "\tcalc 1st order CG_DES Both\n";
	double ret = calcValGradTausAdj( asa->g, asa->x, asa->n );
	return ret;
}

void calcGradASAPiece( asa_objective* asa )
{
	cout << "\tcalc 1st order CG_DES Grad\n";
	calcValGradTausAdj( asa->g, asa->x, asa->n );
}

