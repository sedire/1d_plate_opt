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
	totTimeSteps( CHAR_TIME / DELTA_T + 1 ),
	Km( 0 ),
	dx( 0 ),

	curTime( CHAR_TIME ),
	curTimeStep( totTimeSteps - 1 ),	//CAUTION here
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

	stress_type( 0 ),
	current_type( 0 ),

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
	primSolnDt( 0 )
{
}

AdjSolver::~AdjSolver()
{
	if( primSolnDt != 0 )
	{
		cout << " -- AdjSolver: deleting the dt of the solution array now...\n";
		delete[] primSolnDt;
		cout << " -- AdjSolver: dt of the solution array is deleted!\n";
	}
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
	curTime -= dt;
}

N_PRES AdjSolver::getCurTime()
{
	return curTime;
}

void AdjSolver::loadParamsFromStruct( const SolverPar& loadFrom )
{
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

	stress_type = loadFrom.stress_type;
	current_type = loadFrom.current_type;

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
	for( int i = 0; i < mesh.size(); ++i ){
		mesh[i].setup( eq_num );
	}

	for( int i = 0; i < EQ_NUM; ++i )
	{
		for( int j = 0; j < EQ_NUM; ++j )
		{
			matrA(i, j) = 0.0;
		}
		vectF( i ) = 0.0;
	}

	newmarkA.resize( eq_num, 0.0 );
	newmarkB.resize( eq_num, 0.0 );
}

void AdjSolver::setPrimalSolnData( N_PRES* _primSoln )
{
	if( _primSoln != 0 )
	{
		primSoln = _primSoln;

		//create and fill the array with the approximation of the time derivative of the primal solution
		primSolnDt = new N_PRES[totTimeSteps * Km * eq_num];		//warning here!!!!
		//fill the dt data
		//at the last time step -- backwards finite difference
		for( int y = 0; y < NODES_Y; ++y )
		{
			for( int i = 0; i < EQ_NUM; ++i )
			{
				primSolnDt[( totTimeSteps - 1 ) * Km * eq_num + y * eq_num + i] 
					= ( primSoln[( totTimeSteps - 1 ) * Km * eq_num + y * eq_num + i] - primSoln[( totTimeSteps - 2 ) * Km * eq_num + y * eq_num + i] ) / dt;
			}
		}
		//central finite difference in the middle
		for( int t = totTimeSteps - 2; t >= 1; --t )
		{
			for( int y = 0; y < NODES_Y; ++y )
			{
				for( int i = 0; i < EQ_NUM; ++i )
				{
					primSolnDt[t * Km * eq_num + y * eq_num + i] 
						= ( primSoln[( t + 1 ) * Km * eq_num + y * eq_num + i] - primSoln[( t - 1 ) * Km * eq_num + y * eq_num + i] ) / 2.0 / dt;
				}
			}
		}
		//at the first time step -- forward finite difference
		for( int y = 0; y < NODES_Y; ++y )
		{
			for( int i = 0; i < EQ_NUM; ++i )
			{
				primSolnDt[0 * Km * eq_num + y * eq_num + i] 
					= ( primSoln[1 * Km * eq_num + y * eq_num + i] - primSoln[0 * Km * eq_num + y * eq_num + i] ) / dt;
			}
		}
	}
	else
	{
		cout << "WARNING! pointer to the solution of the primal problem is NULL!\n";
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

void AdjSolver::calcSystemMatrices( int y )
{
	N_PRES Jx = 0.0;
	if( current_type == current_const )
	{
		Jx = J0;
	}
	else if( current_type == current_sin )
	{
		Jx = J0 * sin( (long double)M_PI / tauSin * curTime );
	}
	else if( current_type == current_exp_sin )
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

	matrA( 0, 3 ) = -h * ( 2 * rho + dt * sigma_x * primSoln[indty + 7] * ( primSoln[indty + 7] - 4.0 * beta * dt * primSolnDt[indty + 7] ) ) / ( 2.0 * beta * dt * dt );
	matrA( 0, 4 ) = By1 * h * sigma_x * ( primSoln[indty + 7] - 2.0 * beta * dt * primSolnDt[indty + 7] ) / ( 4.0 * beta * dt );
	matrA( 0, 7 ) = -sigma_x_mu * primSoln[indty + 7] / ( 2.0 * beta * dt ) + sigma_x_mu * primSolnDt[indty + 7];

	matrA( 1, 3 ) = By1 * h * sigma_x * ( primSoln[indty + 7] - 2.0 * beta * dt * primSolnDt[indty + 7] ) / ( 4.0 * beta * dt );
	matrA( 1, 4 ) = -h * ( 8.0 * rho + By1 * By1 * dt * sigma_x ) / ( 8.0 * beta * dt * dt );
	matrA( 1, 7 ) = By1 * sigma_x_mu / ( 4.0 * beta * dt );

	matrA( 2, 1 ) = -1.0;
	matrA( 2, 3 ) = By1 * eps_x_0 * h * ( primSoln[indty + 6] - 2.0 * beta * dt * primSolnDt[indty + 6] ) / ( 4.0 * beta * dt );
	matrA( 2, 4 ) = -eps_x_0 * h * ( 2.0 * beta * dt * primSoln[indty + 6] * primSolnDt[indty + 7] 
		- primSoln[indty + 7] * ( primSoln[indty + 6] - 2.0 * beta * dt * primSolnDt[indty + 6] ) ) / ( 2.0 * beta * dt );
	matrA( 2, 5 ) = h * h * h * ( 2.0 * rho + dt * sigma_x * primSoln[indty + 7] * ( primSoln[indty + 7] - 4.0 * beta * dt * primSolnDt[indty + 7] ) ) 
		/ ( 24.0 * beta * dt * dt );

	matrA( 3, 0 ) = -1.0 / ( B22 * h );
	matrA( 3, 3 ) -eps_x_0 * h * ( -2.0 * beta * dt * primSoln[indty + 6] * primSolnDt[indty + 7] + primSoln[indty + 7] 
		* ( primSoln[indty + 6] - 2.0 * beta * dt * primSolnDt[indty + 6] ) ) / ( 2.0 * B22 * beta * dt );

	matrA( 4, 5 ) = -1.0;

	matrA( 5, 2 ) = 12.0 / ( B22 * h * h * h );
	matrA( 5, 5 ) = -eps_x_0 * ( -2.0 * beta * dt * primSoln[indty + 6] * primSolnDt[indty + 7] + primSoln[indty + 7]
		* ( primSoln[indty + 6] - 2.0 * beta * dt * primSolnDt[indty + 6] ) ) / ( 2.0 * B22 * beta * dt );

	matrA( 6, 3 ) = 0.5 * h * ( By1 * primSolnDt[indty + 2] * eps_x_0 - 2.0 * ( primSolnDt[indty + 3] * eps_x_0 + B22 * sigma_x ) * primSoln[indty + 7] / B22 );
	matrA( 6, 4 ) = 0.5 * By1 * h * sigma_x + primSolnDt[indty + 2] * eps_x_0 * h * primSoln[indty + 7];
	matrA( 6, 5 ) = -primSolnDt[indty + 5] * eps_x_0 * primSoln[indty + 7] / B22;
	matrA( 6, 7 ) = -sigma_x_mu;

	matrA( 7, 3 ) = 0.5 * h * ( -2.0 * Jx + By1 * primSolnDt[indty + 1] * sigma_x - 2.0 * primSolnDt[indty + 3] * eps_x_0 * primSoln[indty + 6] / B22
		- 2.0 * sigma_x * ( 2.0 * primSolnDt[indty + 0] * primSoln[indty + 7] + primSoln[indty + 6] ) );
	matrA( 7, 4 ) = 0.5 * By1 * primSolnDt[indty + 0] * h * sigma_x + primSolnDt[indty + 2] * eps_x_0 * h * primSoln[indty + 6];
	matrA( 7, 5 ) = primSolnDt[indty + 2] * h * h * h * sigma_x * primSoln[indty + 7] / 6.0 - primSolnDt[indty + 5] * eps_x_0 * primSoln[indty + 6] / B22;
	matrA( 7, 6 ) = -1.0 / ( 2.0 * beta * dt );
	matrA( 7, 7 ) = -primSolnDt[indty + 0] * sigma_x_mu;


	vectF( 0 ) = -h * newmarkB[3] * rho + 0.5 * primSoln[indty + 7] 
		* ( -By1 * sigma_x * h * newmarkA[4] + 2.0 * sigma_x_mu * newmarkA[7] + 2.0 * sigma_x * h * newmarkA[3] * primSoln[indty + 7] );
	vectF( 1 ) = 0.25 * ( -4.0 * h * newmarkB[4] * rho + By1 * ( By1 * h * newmarkA[4] * sigma_x - 2.0 * sigma_x_mu * newmarkA[7] ) 
		- 2.0 * By1 * h * newmarkA[3] * sigma_x * primSoln[indty + 7] - 8.0 * primSoln[indty + 1] );
	vectF( 2 ) = 1.0 / 12.0 * ( h * h * h * ( newmarkB[5] * rho - newmarkA[5] * sigma_x * primSoln[indty + 7] * primSoln[indty + 7] )
		- 6.0 * eps_x_0 * h * ( By1 * newmarkA[3] + 2.0 * newmarkA[4] * primSoln[indty + 7] ) * primSoln[indty + 6] );
	vectF( 3 ) = eps_x_0 * h * newmarkA[3] * primSoln[indty + 7] * primSoln[indty + 6] / B22;
	vectF( 5 ) = eps_x_0 * newmarkA[5] * primSoln[indty + 7] * primSoln[indty + 6] / B22;
	vectF( 7 ) = newmarkA[6];
}

N_PRES AdjSolver::doStep()
{
//initial superposition vectors to satisfy the boundary conditions: (WARNING -- the boundary conditions are homogeneous, I'll probably leave the 5th vector as 0)
	N1( 0 ) = 1.0;	N2( 0 ) = 0.0;	N3( 0 ) = 0.0;	N4( 0 ) = 0.0;	N5( 0 ) = 0.0;
	N1( 1 ) = 0.0;	N2( 1 ) = 1.0;	N3( 1 ) = 0.0;	N4( 1 ) = 0.0;	N5( 1 ) = 0.0;
	N1( 2 ) = 0.0;	N2( 2 ) = 0.0;	N3( 2 ) = 0.0;	N4( 2 ) = 0.0;	N5( 2 ) = 0.0;
	N1( 3 ) = 0.0;	N2( 3 ) = 0.0;	N3( 3 ) = 0.0;	N4( 3 ) = 0.0;	N5( 3 ) = 0.0;
	N1( 4 ) = 0.0;	N2( 4 ) = 0.0;	N3( 4 ) = 0.0;	N4( 4 ) = 0.0;	N5( 4 ) = 0.0;
	N1( 5 ) = 0.0;	N2( 5 ) = 0.0;	N3( 5 ) = 1.0;	N4( 5 ) = 0.0;	N5( 5 ) = 0.0;
	N1( 6 ) = 0.0;	N2( 6 ) = 0.0;	N3( 6 ) = 0.0;	N4( 6 ) = 1.0;	N5( 6 ) = 0.0;
	N1( 7 ) = 0.0;	N2( 7 ) = 0.0;	N3( 7 ) = 0.0;	N4( 7 ) = 0.0;	N5( 7 ) = 0.0;

	orthoBuilder->flushO( 0 );
	orthoBuilder->setInitVects( N1, N2, N3, N4, N5 );

	for( int y = 0; y < Km - 1; ++y )
	{
		calcNewmarkAB( y );
		calcSystemMatrices( y );

		rungeKutta->adjCalc( matrA, vectF, dx, 0, &N1 );
		rungeKutta->adjCalc( matrA, vectF, dx, 0, &N2 );
		rungeKutta->adjCalc( matrA, vectF, dx, 0, &N3 );
		rungeKutta->adjCalc( matrA, vectF, dx, 0, &N4 );
		rungeKutta->adjCalc( matrA, vectF, dx, 1, &N5 );

		

		orthoBuilder->flushO( y + 1 );

		orthoBuilder->orthonorm( 1, y, &N1 );
		orthoBuilder->orthonorm( 2, y, &N2 );
		orthoBuilder->orthonorm( 3, y, &N3 );
		orthoBuilder->orthonorm( 4, y, &N4 );
		orthoBuilder->orthonorm( 5, y, &N5 );

		//cout << "\torthog check\n";
		//cout << "\t" << ( N1.transpose() * N2 ) / N1.norm() / N2.norm() << endl;
		//cout << "\t" << ( N1.transpose() * N3 ) / N1.norm() / N3.norm() << endl;
		//cout << "\t" << ( N1.transpose() * N4 ) / N1.norm() / N4.norm() << endl;
		//cout << "\t" << ( N1.transpose() * N5 ) / N1.norm() / N5.norm() << endl;

		//cout << "\t" << ( N2.transpose() * N3 ) / N2.norm() / N3.norm() << endl;
		//cout << "\t" << ( N2.transpose() * N4 ) / N2.norm() / N4.norm() << endl;
		//cout << "\t" << ( N2.transpose() * N5 ) / N2.norm() / N5.norm() << endl;

		//cout << "\t" << ( N3.transpose() * N4 ) / N3.norm() / N4.norm() << endl;
		//cout << "\t" << ( N3.transpose() * N5 ) / N3.norm() / N5.norm() << endl;

		//cout << "\t" << ( N4.transpose() * N5 ) / N4.norm() / N5.norm() << endl;
	}

	cout << "---------------\n";
	cout << N1 << endl;
	cout << "---------------\n";
	cout << N2 << endl;
	cout << "---------------\n";
	cout << N3 << endl;
	cout << "---------------\n";
	cout << N4 << endl;
	cout << "---------------\n";
	cout << N5 << endl;
	cout << "---------------\n";

	return 0;
}