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
	Km( 0 ),
	dx( 0 ),

	curTime( 0 ),
	curTimeStep( 0 ),

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

	eq_num( EQ_NUM )
{
}

void AdjSolver::decreaseTime()
{
	--curTimeStep;
	curTime -= dt;
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

	mesh.resize( 0 );
	mesh.resize( Km );
	for( int i = 0; i < mesh.size(); ++i ){
		mesh[i].setup( eq_num );
	}
}