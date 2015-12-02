#include "Optimizer.h"

void optimizeASA_Taus()
{
	time_t totOptStart = time( 0 );
	cout << "optimizeASA enter\n";

	const N_PRES threshold( 1.e-10 );

	double* x = new double[GRAD_SIZE_FULL];
	for( int i = 0; i < GRAD_SIZE; ++i )
	{
		x[i] = GlobalParams.getCurrentParams1st( i );
	}
	for( int scen = 0; scen < GlobalParams.getNumberOfScenarios(); ++scen )
	{
		for( int i = 0; i < GRAD_SIZE; ++i )
		{
			x[( scen + 1 ) * GRAD_SIZE + i] = GlobalParams.getCurrentParams2nd( scen, i );
		}
	}
	double* lo = new double[GRAD_SIZE_FULL];
	double* hi = new double[GRAD_SIZE_FULL];
	for( int i = 0; i < GRAD_SIZE_FULL; i += 3 )
	{
		lo[i] = -1.0;
		lo[i + 1] = 0.002;
		lo[i + 2] = 0.00001;

		hi[i] = 1.0;
		hi[i + 1] = 1000000000.0;
		hi[i + 2] = 1000000000.0;
	}

	asacg_parm cgParm;
    asa_parm asaParm;
	asa_cg_default( &cgParm );
    asa_default( &asaParm );
    cgParm.PrintParms = TRUE;
    cgParm.PrintLevel = 3;
    asaParm.PrintParms = TRUE;
    asaParm.PrintLevel = 3;

	asa_cg( x, lo, hi, GRAD_SIZE_FULL, NULL, &cgParm, &asaParm, threshold, calcValASA_Taus, calcGradASA_Taus, calc1stOrdOptInfoASA_Taus, 0, 0 );

	cout << "\n\n===============\nASA optimization complete. X is:\n";
	ofstream of( "solution.txt" );
	for( int i = 0; i < GRAD_SIZE_FULL; ++i )
	{
		of << x[i] << endl;
		cout << x[i] << endl;
	}
	of.close();

	cout << " total optimization time : " << time( 0 ) - totOptStart << endl;

	delete[] x;
	delete[] lo;
	delete[] hi;
}

void optimizeASAPiece()
{
	time_t totOptStart = time( 0 );
	cout << "optimizeASA enter\n";

	if( GlobalParams.getCurrentType() != currentPieceLin )
	{
		cout << "ERROR: wrong type of the current\n";
		return;
	}

	const N_PRES threshold( 1.e-10 );

	N_PRES dT = GlobalParams.getTotalTimeT1() / GlobalParams.getNumberOfCurrentParams();
	int firstStageParamNum = (int)floor( GlobalParams.getSwitchTimeT0() / dT );
	int secondStageParamNum = GlobalParams.getNumberOfCurrentParams() - firstStageParamNum;
	int xSize = firstStageParamNum + GlobalParams.getNumberOfScenarios() * secondStageParamNum;
	double* x = new double[xSize];
	cout << " --------------\n\n";
	for( int i = 0; i < firstStageParamNum; ++i )
	{
		x[i] = GlobalParams.getCurrentParams1st( i );
		cout << x[i] << endl;
	}
	for( int scen = 0; scen < GlobalParams.getNumberOfScenarios(); ++scen )
	{
		for( int i = 0; i < secondStageParamNum; ++i )
		{
			x[firstStageParamNum + scen * secondStageParamNum + i] = GlobalParams.getCurrentParams2nd( scen, i );
			cout << x[firstStageParamNum + scen * secondStageParamNum + i] << endl;
		}
	}
	cout << " --------------\n\n";

	double* lo = new double[xSize];
	double* hi = new double[xSize];
	for( int i = 0; i < xSize; ++i )
	{
		lo[i] = -10.0;
		hi[i] = 10.0;
	}

	asacg_parm cgParm;
    asa_parm asaParm;
	asa_cg_default( &cgParm );
    asa_default( &asaParm );
    cgParm.PrintParms = TRUE;
    cgParm.PrintLevel = 3;
    asaParm.PrintParms = TRUE;
    asaParm.PrintLevel = 3;

	asa_cg( x, lo, hi, xSize, NULL, &cgParm, &asaParm, threshold, calcValASA_Taus, calcGradASAPiece, calc1stOrdOptInfoASAPiece, 0, 0 );
	//asa_cg( x, lo, hi, xSize, NULL, &cgParm, &asaParm, threshold, calcValASA_Taus, calcGradASA_Taus, calc1stOrdOptInfoASA_Taus, 0, 0 );

	cout << "\n\n===============\nASA optimization complete. X is:\n";
	ofstream of( "solution.txt" );
	for( int i = 0; i < xSize; ++i )
	{
		of << x[i] << endl;
		cout << x[i] << endl;
	}
	of.close();

	cout << " total optimization time : " << time( 0 ) - totOptStart << endl;

	delete[] x;
	delete[] lo;
	delete[] hi;
}