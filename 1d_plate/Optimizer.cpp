#include "Optimizer.h"


void optimizeASAPiece()
{
	time_t totOptStart = time( 0 );
	cout << "optimizeASA enter\n";

	if( GlobalParams.getCurrentType() != currentPieceLin )
	{
		cout << "ERROR: wrong type of the current\n";
		return;
	}

	const N_PRES threshold( 1.e-6 );

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

	asa_cg( x, lo, hi, GRAD_SIZE_FULL, NULL, &cgParm, &asaParm, threshold, calcValASA_Taus, calcGradASAPiece, calc1stOrdOptInfoASAPiece, 0, 0 );

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