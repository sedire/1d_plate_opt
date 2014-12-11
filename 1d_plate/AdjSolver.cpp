#include "AdjSolver.h"

AdjSolver::AdjSolver() :
	E1( 0 ),			
	E2( 0 ),
	nu21( 0 ),
	rho( 0 ),

	sigma_x( 0 ),
	sigma_x_mu( 0 ),

	h( 0 ),
	a( 0 ),

	dt( 0 ),
	totTimeSteps( ( int )( CHAR_TIME / DELTA_T + 1 ) ),
	Km( 0 ),
	dx( 0 ),

	curTime( CHAR_TIME - DELTA_T ),
	curTimeStep( totTimeSteps - 2 ),	//CAUTION here
	switchTime( SWITCH_TIME ),

	J0( 0 ),
	J0_1( 0 ),
	tauSin( 0 ),
	tauSin_1( 0 ),
	tauExp( 0 ),
	tauExp_1( 0 ),

	p0( 0 ),
	tauP( 0 ),
	rad( 0 ),

	stressType( 0 ),
	currentType( 0 ),

	beta( 0 ),

	B11( 0 ),
	B22( 0 ),
	B12( 0 ),
	By0( 0 ),
	By1( 0 ),

	eps_0( 0 ),
	eps_x( 0 ),
	eps_x_0( 0 ),

	rungeKutta( 0 ),
	orthoBuilder( 0 ),

	eq_num( EQ_NUM ),

	primSoln( 0 ),
	primSolnDt( 0 ),
	adjSoln( 0 )
{
}

AdjSolver::~AdjSolver()
{
	//if( primSolnDt != 0 )
	//{
	//	cout << " -- AdjSolver: deleting the dt of the solution array now...\n";
	//	delete[] primSolnDt;
	//	cout << " -- AdjSolver: dt of the solution array is deleted!\n";
	//}
	if( rungeKutta != 0 )
	{
		delete rungeKutta;
	}
	if( orthoBuilder != 0 )
	{
		delete orthoBuilder;
	}
}

void AdjSolver::decreaseTime()
{
	--curTimeStep;
	if( curTimeStep == 0 )
	{
		curTime = 0.0;
	}
	else
	{
		curTime -= dt;
	}
}

N_PRES AdjSolver::getCurTime()
{
	return curTime;
}

int AdjSolver::getCurTimeStep()
{
	return curTimeStep;
}

void AdjSolver::loadParamsFromStruct( const SolverPar& loadFrom )
{
	ofstream of1( "adj_test_sol.txt" );
	of1.close();

	E1 = loadFrom.E1;				
	E2 = loadFrom.E2;				
	nu21 = loadFrom.nu21;			
	rho = loadFrom.rho;				

	sigma_x = loadFrom.sigma_x;			
	sigma_x_mu = loadFrom.sigma_x_mu;

	h = loadFrom.h;
	a = loadFrom.a;

	dt = loadFrom.dt;
	Km = loadFrom.Km;
	dx = loadFrom.dx;

	J0 = loadFrom.J0;
	J0_1 = loadFrom.J0_1;
	tauSin = loadFrom.tauSin;
	tauSin_1 = loadFrom.tauSin_1;
	tauExp = loadFrom.tauExp;
	tauExp_1 = loadFrom.tauExp_1;

	p0 = loadFrom.p0;
	tauP = loadFrom.tauP;
	rad = loadFrom.rad;

	stressType = loadFrom.stressType;
	currentType = loadFrom.currentType;

	beta = loadFrom.beta;

	B11 = loadFrom.B11;
	B22 = loadFrom.B22;
	B12 = loadFrom.B12;
	By0 = loadFrom.By0;
	By1 = loadFrom.By1;                                      								

	eps_0 = loadFrom.eps_0;
	eps_x = loadFrom.eps_x;
	eps_x_0 = loadFrom.eps_x_0;

	if( rungeKutta != 0 )
	{
		delete rungeKutta;
	}
	rungeKutta = new RungeKutta<N_PRES>( eq_num );
	if( orthoBuilder != 0 )
	{
		delete orthoBuilder;
	}
	orthoBuilder = new OrthoBuilderGSh<N_PRES>();
	orthoBuilder->setParams( Km );

//resetting all data containers to 0s
	mesh.resize( 0 );
	mesh.resize( Km );
	for( int i = 0; i < mesh.size(); ++i )
	{
		mesh[i].setup( eq_num );
	}

	for( int i = 0; i < EQ_NUM; ++i )
	{
		for( int j = 0; j < EQ_NUM; ++j )
		{
			matrA( i, j ) = 0.0;
			matrA1( i, j ) = 0.0;
		}
		vectF( i ) = 0.0;
		vectF1( i ) = 0.0;

		N1( i ) = 0.0l;
		N2( i ) = 0.0l;
		N3( i ) = 0.0l;
		N4( i ) = 0.0l;
		N5( i ) = 0.0l;

		N1orthog( i ) = 0.0l;
		N2orthog( i ) = 0.0l;
		N3orthog( i ) = 0.0l;
		N4orthog( i ) = 0.0l;
		N5orthog( i ) = 0.0l;
	}

	newmarkA.resize( eq_num, 0.0 );
	newmarkB.resize( eq_num, 0.0 );
}

void AdjSolver::setPrimalSolnData( N_PRES* _primSoln )
{
	if( _primSoln != 0 )
	{
		primSoln = _primSoln;

		////create and fill the array with the approximation of the time derivative of the primal solution
		//primSolnDt = new N_PRES[totTimeSteps * Km * eq_num];		//warning here!!!!
		////fill the dt data
		////at the last time step -- backwards finite difference
		//for( int y = 0; y < NODES_Y; ++y )
		//{
		//	for( int i = 0; i < EQ_NUM; ++i )
		//	{
		//		primSolnDt[( totTimeSteps - 1 ) * Km * eq_num + y * eq_num + i] 
		//			= ( primSoln[( totTimeSteps - 1 ) * Km * eq_num + y * eq_num + i] - primSoln[( totTimeSteps - 2 ) * Km * eq_num + y * eq_num + i] ) / dt;
		//	}
		//}
		////central finite difference in the middle
		//for( int t = totTimeSteps - 2; t >= 1; --t )
		//{
		//	for( int y = 0; y < NODES_Y; ++y )
		//	{
		//		for( int i = 0; i < EQ_NUM; ++i )
		//		{
		//			primSolnDt[t * Km * eq_num + y * eq_num + i] 
		//				= ( primSoln[( t + 1 ) * Km * eq_num + y * eq_num + i] - primSoln[( t - 1 ) * Km * eq_num + y * eq_num + i] ) / 2.0 / dt;
		//		}
		//	}
		//}
		////at the first time step -- forward finite difference
		//for( int y = 0; y < NODES_Y; ++y )
		//{
		//	for( int i = 0; i < EQ_NUM; ++i )
		//	{
		//		primSolnDt[0 * Km * eq_num + y * eq_num + i] 
		//			= 0.0/*( primSoln[1 * Km * eq_num + y * eq_num + i] - primSoln[0 * Km * eq_num + y * eq_num + i] ) / dt*/;
		//	}
		//}
	}
	else
	{
		cout << "WARNING! pointer to the solution of the primal problem is NULL!\n";
	}
}

void AdjSolver::setPrimalDtData( N_PRES* _primSolnDt )
{
	if( _primSolnDt != 0 )
	{
		primSolnDt = _primSolnDt;
	}
	else
	{
		cout << "WARNING! pointer to the Dt of the solution of the primal problem is NULL!\n";
	}
}


void AdjSolver::setAdjointSolnData( N_PRES* _adjSoln )
{
	if( _adjSoln != 0 )
	{
		adjSoln = _adjSoln;
		for( int i = 0; i < ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM; ++i )
		{
			adjSoln[i] = 0.0;
		}
	}
	else
	{
		cout << "WARNING! pointer to the solution of the ajoint problem is NULL!\n";
	}
}

void AdjSolver::calcNewmarkAB( int y )
{
	for( int i = 0; i < eq_num; ++i )
	{
		newmarkA[i] = 1.0 / ( 2.0 * beta * dt ) * mesh[y].N1[i] - 2.0 * ( 0.5 - beta ) / ( 2.0 * beta ) * mesh[y].d1N1[i] 
			- ( 2.0 * beta - 0.5 ) * dt / ( 2.0 * beta ) * mesh[y].d2N1[i];
		newmarkB[i] = -1.0 / ( beta * dt * dt ) * mesh[y].N1[i] + 1.0 / ( beta * dt ) * mesh[y].d1N1[i]
			+ ( 1.0 - 1.0 / ( 2.0 * beta ) ) * mesh[y].d2N1[i];
	}
}

void AdjSolver::calcSystemMatrices( int y, Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor>* A, Matrix<PL_NUM, EQ_NUM, 1>* f )
{
	N_PRES Jx = 0.0;
	if( currentType == current_const )
	{
		Jx = J0;
	}
	else if( currentType == current_sin )
	{
		Jx = J0 * sin( (long double)M_PI / tauSin * curTime );
	}
	else if( currentType == current_exp_sin )
	{
		if( curTime <= switchTime )
		{
			Jx = J0 * exp( -curTime / tauExp ) * sin( (long double)M_PI / tauSin * curTime );
		}
		else
		{
			Jx = J0_1 * exp( -curTime / tauExp_1 ) * sin( (long double)M_PI / tauSin_1 * curTime );
		}
	}

	int indty = curTimeStep * Km * eq_num + y * eq_num;

	(*A)( 0, 3 ) = -h * ( 2 * rho + dt * sigma_x * primSoln[indty + 7] * ( primSoln[indty + 7] - 4.0 * beta * dt * primSolnDt[indty + 7] ) ) / ( 2.0 * beta * dt * dt );
	(*A)( 0, 4 ) = By1 * h * sigma_x * ( primSoln[indty + 7] - 2.0 * beta * dt * primSolnDt[indty + 7] ) / ( 4.0 * beta * dt );
	(*A)( 0, 7 ) = -sigma_x_mu * primSoln[indty + 7] / ( 2.0 * beta * dt ) + sigma_x_mu * primSolnDt[indty + 7];

	(*A)( 1, 3 ) = By1 * h * sigma_x * ( primSoln[indty + 7] - 2.0 * beta * dt * primSolnDt[indty + 7] ) / ( 4.0 * beta * dt );
	(*A)( 1, 4 ) = -h * ( 8.0 * rho + By1 * By1 * dt * sigma_x ) / ( 8.0 * beta * dt * dt );
	(*A)( 1, 7 ) = By1 * sigma_x_mu / ( 4.0 * beta * dt );

	(*A)( 2, 1 ) = -1.0;
	(*A)( 2, 3 ) = By1 * eps_x_0 * h * ( primSoln[indty + 6] - 2.0 * beta * dt * primSolnDt[indty + 6] ) / ( 4.0 * beta * dt );
	(*A)( 2, 4 ) = 0.5 * eps_x_0 * h * ( -2.0 * primSoln[indty + 6] * primSolnDt[indty + 7] + primSoln[indty + 7] * ( primSoln[indty + 6] / ( beta * dt ) 
		-2.0 * primSolnDt[indty + 6] ) );
	(*A)( 2, 5 ) = h * h * h * ( 2.0 * rho + dt * sigma_x * primSoln[indty + 7] * ( primSoln[indty + 7] - 4.0 * beta * dt * primSolnDt[indty + 7] ) ) 
		/ ( 24.0 * beta * dt * dt );

	(*A)( 3, 0 ) = -1.0 / ( B22 * h );
	(*A)( 3, 3 ) = eps_x_0 * h * ( 2.0 * beta * dt * primSoln[indty + 6] * primSolnDt[indty + 7] - primSoln[indty + 7] 
		* ( primSoln[indty + 6] - 2.0 * beta * dt * primSolnDt[indty + 6] ) ) / ( 2.0 * B22 * beta * dt );

	(*A)( 4, 5 ) = -1.0;

	(*A)( 5, 2 ) = 12.0 / ( B22 * h * h * h );
	(*A)( 5, 5 ) = eps_x_0 * ( 2.0 * beta * dt * primSoln[indty + 6] * primSolnDt[indty + 7] - primSoln[indty + 7]
		* ( primSoln[indty + 6] - 2.0 * beta * dt * primSolnDt[indty + 6] ) ) / ( 2.0 * B22 * beta * dt );

	(*A)( 6, 3 ) = 0.5 * By1 * primSolnDt[indty + 2] * eps_x_0 * h 
		- h * ( primSolnDt[indty + 3] * eps_x_0 + B22 * sigma_x ) * primSoln[indty + 7] / B22;
	(*A)( 6, 4 ) = 0.5 * By1 * h * sigma_x + primSolnDt[indty + 2] * eps_x_0 * h * primSoln[indty + 7];
	(*A)( 6, 5 ) = -primSolnDt[indty + 5] * eps_x_0 * primSoln[indty + 7] / B22;
	(*A)( 6, 7 ) = -sigma_x_mu;

	(*A)( 7, 3 ) = 0.5 * h * ( -2.0 * Jx + By1 * primSolnDt[indty + 1] * sigma_x - 4.0 * primSolnDt[indty + 0] * sigma_x * primSoln[indty + 7]
		- 2.0 * ( primSolnDt[indty + 3] * eps_x_0 + B22 * sigma_x ) * primSoln[indty + 6] / B22 );
	(*A)( 7, 4 ) = 0.5 * By1 * primSolnDt[indty + 0] * h * sigma_x + primSolnDt[indty + 2] * eps_x_0 * h * primSoln[indty + 6];
	(*A)( 7, 5 ) = primSolnDt[indty + 2] * h * h * h * sigma_x * primSoln[indty + 7] / 6.0 - primSolnDt[indty + 5] * eps_x_0 * primSoln[indty + 6] / B22;
	(*A)( 7, 6 ) = -1.0 / ( 2.0 * beta * dt );
	(*A)( 7, 7 ) = -primSolnDt[indty + 0] * sigma_x_mu;


	(*f)( 0 ) = -h * newmarkB[3] * rho + primSoln[indty + 7] 
		* ( -0.5 * By1 * sigma_x * h * newmarkA[4] + sigma_x_mu * newmarkA[7] + sigma_x * h * newmarkA[3] * primSoln[indty + 7] );
	(*f)( 1 ) = 0.25 * ( -4.0 * h * newmarkB[4] * rho + By1 * ( By1 * h * newmarkA[4] * sigma_x - 2.0 * sigma_x_mu * newmarkA[7] ) 
		- 2.0 * By1 * h * newmarkA[3] * sigma_x * primSoln[indty + 7] - 8.0 * primSoln[indty + 1] );
	(*f)( 2 ) = 1.0 / 12.0 * ( h * h * h * ( newmarkB[5] * rho - newmarkA[5] * sigma_x * primSoln[indty + 7] * primSoln[indty + 7] )
		- 6.0 * eps_x_0 * h * ( By1 * newmarkA[3] + 2.0 * newmarkA[4] * primSoln[indty + 7] ) * primSoln[indty + 6] );
	(*f)( 3 ) = eps_x_0 * h * newmarkA[3] * primSoln[indty + 7] * primSoln[indty + 6] / B22;
	(*f)( 5 ) = eps_x_0 * newmarkA[5] * primSoln[indty + 7] * primSoln[indty + 6] / B22;
	(*f)( 7 ) = newmarkA[6];
}

//void AdjSolver::calcSystemMatrices( int y )	//This one is scaled. I did this because in some of my computations one of the components
												//of the solution of the adjoint behaved weirdly. I wanted to scale it up to see what happens
												//(because I thought that is caused by a numerical error from numbers being too small)
												//but in fact it doesn't seem it influences it in any way.
//{
//	N_PRES coef = 1;
//	N_PRES Jx = 0.0;
//	if( currentType == current_const )
//	{
//		Jx = J0;
//	}
//	else if( currentType == current_sin )
//	{
//		Jx = J0 * sin( (long double)M_PI / tauSin * curTime );
//	}
//	else if( currentType == current_exp_sin )
//	{
//		if( curTime <= switchTime )
//		{
//			Jx = J0 * exp( -curTime / tauExp ) * sin( (long double)M_PI / tauSin * curTime );
//		}
//		else
//		{
//			Jx = J0_1 * exp( -curTime / tauExp_1 ) * sin( (long double)M_PI / tauSin_1 * curTime );
//		}
//	}
//
//	int indty = curTimeStep * Km * eq_num + y * eq_num;
//
//	matrA( 0, 3 ) = ( -h * ( 2 * rho + dt * sigma_x * primSoln[indty + 7] * ( primSoln[indty + 7] - 4.0 * beta * dt * primSolnDt[indty + 7] ) ) / ( 2.0 * beta * dt * dt ) ) / coef;
//	matrA( 0, 4 ) = By1 * h * sigma_x * ( primSoln[indty + 7] - 2.0 * beta * dt * primSolnDt[indty + 7] ) / ( 4.0 * beta * dt );
//	matrA( 0, 7 ) = -sigma_x_mu * primSoln[indty + 7] / ( 2.0 * beta * dt ) + sigma_x_mu * primSolnDt[indty + 7];
//
//	matrA( 1, 3 ) = ( By1 * h * sigma_x * ( primSoln[indty + 7] - 2.0 * beta * dt * primSolnDt[indty + 7] ) / ( 4.0 * beta * dt ) ) / coef;
//	matrA( 1, 4 ) = -h * ( 8.0 * rho + By1 * By1 * dt * sigma_x ) / ( 8.0 * beta * dt * dt );
//	matrA( 1, 7 ) = By1 * sigma_x_mu / ( 4.0 * beta * dt );
//
//	matrA( 2, 1 ) = -1.0;
//	matrA( 2, 3 ) = ( By1 * eps_x_0 * h * ( primSoln[indty + 6] - 2.0 * beta * dt * primSolnDt[indty + 6] ) / ( 4.0 * beta * dt ) ) / coef;
//	matrA( 2, 4 ) = -eps_x_0 * h * ( 2.0 * beta * dt * primSoln[indty + 6] * primSolnDt[indty + 7] 
//		- primSoln[indty + 7] * ( primSoln[indty + 6] - 2.0 * beta * dt * primSolnDt[indty + 6] ) ) / ( 2.0 * beta * dt );
//	matrA( 2, 5 ) = h * h * h * ( 2.0 * rho + dt * sigma_x * primSoln[indty + 7] * ( primSoln[indty + 7] - 4.0 * beta * dt * primSolnDt[indty + 7] ) ) 
//		/ ( 24.0 * beta * dt * dt );
//
//	matrA( 3, 0 ) = -1.0 / ( B22 * h ) * coef;
//	matrA( 3, 3 ) = -eps_x_0 * h * ( -2.0 * beta * dt * primSoln[indty + 6] * primSolnDt[indty + 7] + primSoln[indty + 7] 
//		* ( primSoln[indty + 6] - 2.0 * beta * dt * primSolnDt[indty + 6] ) ) / ( 2.0 * B22 * beta * dt );
//
//	matrA( 4, 5 ) = -1.0;
//
//	matrA( 5, 2 ) = 12.0 / ( B22 * h * h * h );
//	matrA( 5, 5 ) = -eps_x_0 * ( -2.0 * beta * dt * primSoln[indty + 6] * primSolnDt[indty + 7] + primSoln[indty + 7]
//		* ( primSoln[indty + 6] - 2.0 * beta * dt * primSolnDt[indty + 6] ) ) / ( 2.0 * B22 * beta * dt );
//
//	matrA( 6, 3 ) = ( 0.5 * h * ( By1 * primSolnDt[indty + 2] * eps_x_0 - 2.0 * ( primSolnDt[indty + 3] * eps_x_0 + B22 * sigma_x ) * primSoln[indty + 7] / B22 ) ) / coef;
//	matrA( 6, 4 ) = 0.5 * By1 * h * sigma_x + primSolnDt[indty + 2] * eps_x_0 * h * primSoln[indty + 7];
//	matrA( 6, 5 ) = -primSolnDt[indty + 5] * eps_x_0 * primSoln[indty + 7] / B22;
//	matrA( 6, 7 ) = -sigma_x_mu;
//
//	matrA( 7, 3 ) = ( 0.5 * h * ( -2.0 * Jx + By1 * primSolnDt[indty + 1] * sigma_x - 2.0 * primSolnDt[indty + 3] * eps_x_0 * primSoln[indty + 6] / B22
//		- 2.0 * sigma_x * ( 2.0 * primSolnDt[indty + 0] * primSoln[indty + 7] + primSoln[indty + 6] ) ) ) / coef;
//	matrA( 7, 4 ) = 0.5 * By1 * primSolnDt[indty + 0] * h * sigma_x + primSolnDt[indty + 2] * eps_x_0 * h * primSoln[indty + 6];
//	matrA( 7, 5 ) = primSolnDt[indty + 2] * h * h * h * sigma_x * primSoln[indty + 7] / 6.0 - primSolnDt[indty + 5] * eps_x_0 * primSoln[indty + 6] / B22;
//	matrA( 7, 6 ) = -1.0 / ( 2.0 * beta * dt );
//	matrA( 7, 7 ) = -primSolnDt[indty + 0] * sigma_x_mu;
//
//
//	vectF( 0 ) = ( -2.0 * h * newmarkB[3] * rho + primSoln[indty + 7] 
//		* ( -By1 * coef * sigma_x * h * newmarkA[4] + 2.0 * coef * sigma_x_mu * newmarkA[7] + 2.0 * sigma_x * h * newmarkA[3] * primSoln[indty + 7] ) ) / ( 2.0 * coef );
//	vectF( 1 ) = 0.25 * ( -4.0 * h * newmarkB[4] * rho + By1 * ( By1 * h * newmarkA[4] * sigma_x - 2.0 * sigma_x_mu * newmarkA[7] ) 
//		- 2.0 * By1 * h * newmarkA[3] * sigma_x * primSoln[indty + 7] / coef - 8.0 * primSoln[indty + 1] );
//	vectF( 2 ) = 1.0 / ( 12.0 * coef ) * ( coef * h * h * h * ( newmarkB[5] * rho - newmarkA[5] * sigma_x * primSoln[indty + 7] * primSoln[indty + 7] )
//		- 6.0 * eps_x_0 * h * ( By1 * newmarkA[3] + 2.0 * coef * newmarkA[4] * primSoln[indty + 7] ) * primSoln[indty + 6] );
//	vectF( 3 ) = eps_x_0 * h * newmarkA[3] * primSoln[indty + 7] * primSoln[indty + 6] / B22;
//	vectF( 5 ) = eps_x_0 * newmarkA[5] * primSoln[indty + 7] * primSoln[indty + 6] / B22;
//	vectF( 7 ) = newmarkA[6];
//}

N_PRES AdjSolver::doStep()
{
//initial superposition vectors to satisfy the boundary conditions: (WARNING -- the boundary conditions are homogeneous, I'll probably leave the 5th vector as 0)
	N1( 0 ) = 1.0;	N2( 0 ) = 0.0;	N3( 0 ) = 0.0;	N4( 0 ) = 0.0;	N5( 0 ) = 0.0;
	N1( 1 ) = 0.0;	N2( 1 ) = 1.0;	N3( 1 ) = 0.0;	N4( 1 ) = 0.0;	N5( 1 ) = 0.0;
	N1( 2 ) = 0.0;	N2( 2 ) = 0.0;	N3( 2 ) = 0.0;	N4( 2 ) = 0.0;	N5( 2 ) = 0.0;
	N1( 3 ) = 0.0;	N2( 3 ) = 0.0;	N3( 3 ) = 0.0;	N4( 3 ) = 0.0;	N5( 3 ) = 0.0;
	N1( 4 ) = 0.0;	N2( 4 ) = 0.0;	N3( 4 ) = 0.0;	N4( 4 ) = 0.0;	N5( 4 ) = 0.0;
	N1( 5 ) = 0.0;	N2( 5 ) = 0.0;	N3( 5 ) = 1.0;	N4( 5 ) = 0.0;	N5( 5 ) = 0.5;
	N1( 6 ) = 0.0;	N2( 6 ) = 0.0;	N3( 6 ) = 0.0;	N4( 6 ) = 1.0;	N5( 6 ) = 0.5;
	N1( 7 ) = 0.0;	N2( 7 ) = 0.0;	N3( 7 ) = 0.0;	N4( 7 ) = 0.0;	N5( 7 ) = 0.0;

	orthoBuilder->flushO( 0 );
	orthoBuilder->setInitVects( N1, N2, N3, N4, N5 );

	int active = 1;		//how many nodes in a row we have, on which orthonormalization was not performed. 
						//We need this to know whem we can switch to ABM method

	calcNewmarkAB( 0 );
	calcSystemMatrices( 0, &matrA1, &vectF1 );

	for( int y = 0; y < Km - 1; ++y )
	{
		matrA = matrA1;
		vectF = vectF1;

		calcNewmarkAB( y + 1 );
		calcSystemMatrices( y + 1, &matrA1, &vectF1 );

		rungeKutta->adjCalc( matrA, matrA1, vectF, vectF1, dx, 0, &N1 );
		rungeKutta->adjCalc( matrA, matrA1, vectF, vectF1, dx, 0, &N2 );
		rungeKutta->adjCalc( matrA, matrA1, vectF, vectF1, dx, 0, &N3 );
		rungeKutta->adjCalc( matrA, matrA1, vectF, vectF1, dx, 0, &N4 );
		rungeKutta->adjCalc( matrA, matrA1, vectF, vectF1, dx, 1, &N5 );

		orthoBuilder->flushO( y + 1 );

		orthoBuilder->orthonorm( 1, y, &N1 );
		orthoBuilder->orthonorm( 2, y, &N2 );
		orthoBuilder->orthonorm( 3, y, &N3 );
		orthoBuilder->orthonorm( 4, y, &N4 );
		orthoBuilder->orthonorm( 5, y, &N5 );
	}

	orthoBuilder->buildSolutionAdj( &mesh );

	for( int x = 0; x < Km; ++x )
	{
		for( int i = 0; i < eq_num; ++i )
		{
			N_PRES d2N = ( mesh[x].N[i] - mesh[x].N1[i] ) / beta / dt / dt + mesh[x].d1N1[i] / beta / dt
						+ ( 1.0 - 1.0 / ( 2.0 * beta ) ) * mesh[x].d2N1[i];
			N_PRES d1N = -0.5 * dt * d2N + mesh[x].d1N1[i] - 0.5 * dt * mesh[x].d2N1[i];
			mesh[x].N1[i] = mesh[x].N[i];
			mesh[x].d1N1[i] = d1N;
			mesh[x].d2N1[i] = d2N;
		}
	}

	//cout << "---------------\n";
	//cout << N1 << endl;
	//cout << "---------------\n";
	//cout << N2 << endl;
	//cout << "---------------\n";
	//cout << N3 << endl;
	//cout << "---------------\n";
	//cout << N4 << endl;
	//cout << "---------------\n";
	//cout << N5 << endl;
	//cout << "---------------\n";

	if( adjSoln != 0 )
	{
		for( int y = 0; y < Km; ++y )
		{
			for( int i = 0; i < eq_num; ++i )
			{
				adjSoln[curTimeStep * ( Km * eq_num ) + y * eq_num + i] = mesh[y].N[i];
			}
		}
	}

	return 0;
}

void AdjSolver::scalePrevTimeStepSoln()
{
	for( int x = 0; x < Km; ++x )
	{
		for( int i = 0; i < eq_num; ++i )
		{
			mesh[x].N1[i] *= SWITCH_TIME / ( CHAR_TIME - SWITCH_TIME );
			mesh[x].d1N1[i] *= SWITCH_TIME / ( CHAR_TIME - SWITCH_TIME );
			mesh[x].d2N1[i] *= SWITCH_TIME / ( CHAR_TIME - SWITCH_TIME );
		}
	}
}

N_PRES AdjSolver::calcJDeriv()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			for( int y = 0; y < Km; ++y )
			{
				sum += exp( -t * dt / tauExp ) * sin( M_PI * t * dt / tauSin ) 
					* ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7] 
					- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 ) * dx * dt; 
			}
		}
	}

	return sum * J0_SCALE;
}

N_PRES AdjSolver::calcJ0Deriv()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			if( t * dt <= SWITCH_TIME )
			{
				for( int y = 0; y < Km; ++y )
				{
					sum += exp( -t * dt / tauExp ) * sin( M_PI * t * dt / tauSin ) 
						* ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7] 
						- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 ) * dx * dt; 
				}
			}
		}
	}

	return sum * J0_SCALE;
}

N_PRES AdjSolver::calcJ0DerivS()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			if( t * dt < SWITCH_TIME )
			{
				for( int y = 0; y < Km; ++y )
				{
					sum += exp( -t * dt / tauExp ) * sin( M_PI * t * dt / tauSin ) 
						* ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7] 
						- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 ) * dx * dt / SWITCH_TIME; 
				}
			}
			else if( t * dt == SWITCH_TIME )
			{
				N_PRES coef = SWITCH_TIME / ( CHAR_TIME - SWITCH_TIME );
				for( int y = 0; y < Km; ++y )
				{
					sum += exp( -t * dt / tauExp ) * sin( M_PI * t * dt / tauSin ) 
						* ( coef * adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7] 
						- coef * adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 ) * dx * dt / SWITCH_TIME; 
				}
			}
		}
	}

	return sum * J0_SCALE;
}

N_PRES AdjSolver::calcJ1Deriv()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			if( t * dt > SWITCH_TIME )
			{
				for( int y = 0; y < Km; ++y )
				{
					sum += exp( -t * dt / tauExp_1 ) * sin( M_PI * t * dt / tauSin_1 ) 
						* ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7] 
						- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 ) * dx * dt; 
				}
			}
		}
	}

	return sum * J0_SCALE;
}

N_PRES AdjSolver::calcJ1DerivS()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			if( t * dt > SWITCH_TIME )
			{
				for( int y = 0; y < Km; ++y )
				{
					sum += exp( -t * dt / tauExp_1 ) * sin( M_PI * t * dt / tauSin_1 ) 
						* ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7] 
						- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 ) * dx * dt / ( CHAR_TIME - SWITCH_TIME ); 
				}
			}
		}
	}

	return sum * J0_SCALE;
}

N_PRES AdjSolver::calcTauSinDeriv()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			N_PRES innerSum = 0.0;
			for( int y = 0; y < Km; ++y )
			{
				innerSum += ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7]
					- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 );
			}
			sum += -innerSum * exp( -t * dt / tauExp ) * cos( t * dt * M_PI / tauSin ) * t * dt;
		}
		sum *= J0 * dx * dt * M_PI / ( tauSin * tauSin );
	}

	return sum;
}

N_PRES AdjSolver::calcTauSin0Deriv()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			if( t * dt <= SWITCH_TIME )
			{
				N_PRES innerSum = 0.0;
				for( int y = 0; y < Km; ++y )
				{
					innerSum += ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7]
						- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 );
				}
				sum += -innerSum * exp( -t * dt / tauExp ) * cos( t * dt * M_PI / tauSin ) * t * dt;
			}
		}
		sum *= J0 * dx * dt * M_PI / ( tauSin * tauSin );
	}

	return sum;
}

N_PRES AdjSolver::calcTauSin0DerivS()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			if( t * dt < SWITCH_TIME )
			{
				N_PRES innerSum = 0.0;
				for( int y = 0; y < Km; ++y )
				{
					innerSum += ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7]
						- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 );
				}
				sum += -innerSum * exp( -t * dt / tauExp ) * cos( t * dt * M_PI / tauSin ) * t * dt;
			}
			else if( t * dt == SWITCH_TIME )
			{
				N_PRES coef = SWITCH_TIME / ( CHAR_TIME - SWITCH_TIME );
				N_PRES innerSum = 0.0;
				for( int y = 0; y < Km; ++y )
				{
					innerSum += ( coef * adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7]
						- coef * adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 );
				}
				sum += -innerSum * exp( -t * dt / tauExp ) * cos( t * dt * M_PI / tauSin ) * t * dt;
			}
		}
		sum *= J0 * dx * dt * M_PI / ( tauSin * tauSin ) / SWITCH_TIME;
	}

	return sum;
}

N_PRES AdjSolver::calcTauSin1Deriv()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			if( t * dt > SWITCH_TIME )
			{
				N_PRES innerSum = 0.0;
				for( int y = 0; y < Km; ++y )
				{
					innerSum += ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7]
						- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 );
				}
				sum += -innerSum * exp( -t * dt / tauExp_1 ) * cos( t * dt * M_PI / tauSin_1 ) * t * dt;
			}
		}
		sum *= J0_1 * dx * dt * M_PI / ( tauSin_1 * tauSin_1 );
	}

	return sum;
}

N_PRES AdjSolver::calcTauSin1DerivS()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			if( t * dt > SWITCH_TIME )
			{
				N_PRES innerSum = 0.0;
				for( int y = 0; y < Km; ++y )
				{
					innerSum += ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7]
						- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 );
				}
				sum += -innerSum * exp( -t * dt / tauExp_1 ) * cos( t * dt * M_PI / tauSin_1 ) * t * dt;
			}
		}
		sum *= J0_1 * dx * dt * M_PI / ( tauSin_1 * tauSin_1 ) / ( CHAR_TIME - SWITCH_TIME );
	}

	return sum;
}

N_PRES AdjSolver::calcTauExpDeriv()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			for( int y = 0; y < Km; ++y )
			{
				sum += J0 * exp( -t * dt / tauExp ) * sin( M_PI * t * dt / tauSin ) * t * dt / ( tauExp * tauExp )
					* ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7] 
					- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 ) * dx * dt;
			}
		}
	}

	return sum;
}

N_PRES AdjSolver::calcTauExp0Deriv()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			if( t * dt <= SWITCH_TIME )
			{
				for( int y = 0; y < Km; ++y )
				{
					sum += J0 * exp( -t * dt / tauExp ) * sin( M_PI * t * dt / tauSin ) * t * dt / ( tauExp * tauExp )
						* ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7] 
						- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 ) * dx * dt;
				}
			}
		}
	}

	return sum;
}

N_PRES AdjSolver::calcTauExp0DerivS()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			if( t * dt < SWITCH_TIME )
			{
				for( int y = 0; y < Km; ++y )
				{
					sum += J0 * exp( -t * dt / tauExp ) * sin( M_PI * t * dt / tauSin ) * t * dt / ( tauExp * tauExp )
						* ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7] 
						- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 ) * dx * dt / SWITCH_TIME;
				}
			}
			else if( t * dt == SWITCH_TIME )
			{
				N_PRES coef = SWITCH_TIME / ( CHAR_TIME - SWITCH_TIME );
				for( int y = 0; y < Km; ++y )
				{
					sum += J0 * exp( -t * dt / tauExp ) * sin( M_PI * t * dt / tauSin ) * t * dt / ( tauExp * tauExp )
						* ( coef * adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7] 
						- coef * adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 ) * dx * dt / SWITCH_TIME;
				}
			}
		}
	}

	return sum;
}

N_PRES AdjSolver::calcTauExp1Deriv()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			if( t * dt > SWITCH_TIME )
			{
				for( int y = 0; y < Km; ++y )
				{
					sum += J0_1 * exp( -t * dt / tauExp_1 ) * sin( M_PI * t * dt / tauSin_1 ) * t * dt / ( tauExp_1 * tauExp_1 )
						* ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7] 
						- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 ) * dx * dt;
				}
			}
		}
	}

	return sum;
}

N_PRES AdjSolver::calcTauExp1DerivS()
{
	N_PRES sum = 0.0;
	
	if( adjSoln != 0 )
	{
		for( int t = 0; t < totTimeSteps; ++t )
		{
			if( t * dt > SWITCH_TIME )
			{
				for( int y = 0; y < Km; ++y )
				{
					sum += J0_1 * exp( -t * dt / tauExp_1 ) * sin( M_PI * t * dt / tauSin_1 ) * t * dt / ( tauExp_1 * tauExp_1 )
						* ( adjSoln[t * ( Km * eq_num ) + y * eq_num + 3] * h * primSoln[t * ( Km * eq_num ) + y * eq_num + 7] 
						- adjSoln[t * ( Km * eq_num ) + y * eq_num + 4] * 0.5 * h * By1 ) * dx * dt / ( CHAR_TIME - SWITCH_TIME );
				}
			}
		}
	}

	return sum;
}

void AdjSolver::dumpSol( int fNum )
{
	N_PRES t = curTime;

	stringstream ss;
	if( fNum >= 0 )
	{
		ss << "adj_test_sol_" << fNum << ".txt";
	}
	else
	{
		ss << "adj_test_sol.txt";
	}
	ofstream of1( ss.str(), ofstream::app );
	of1 << curTime;
	for( int i = 0; i < eq_num; ++i )
	{
		of1 << " ; " << mesh[ ( Km - 1 ) / 2 ].N[i];
	}
	of1 << endl;
	of1.close();
}