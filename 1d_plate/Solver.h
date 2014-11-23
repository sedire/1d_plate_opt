#ifndef _PLATE_1D_SOLVER_
#define _PLATE_1D_SOLVER_ 1

#include <time.h>
#include "plate_var_types.h"
#include "RungeKutta.h"
#include "OrthoBuilder.h"
#include "VarVect.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <omp.h>
#include <Eigen\Eigen>
#include <complex>
#include "SolverPar.h"

using std::cout;
using std::vector;
using std::ofstream;
using std::string;
using std::stringstream;
using std::complex;
using namespace Eigen;

template<class PL_NUM>
class Solver
{
public:
	Solver();
	~Solver();

	int getMaxNewtonIterReached();

	N_PRES* resArr;		//array to keep all the results obtained
	N_PRES* resArrDt;

	void setTask( PL_NUM _J0, PL_NUM _tauSin, PL_NUM _tauExp,
				PL_NUM _J0_1, PL_NUM _tauSin_1, PL_NUM _tauExp_1,
				PL_NUM _By0, int _stressType, PL_NUM _p0, PL_NUM _tauP );
	void setMechLoad( int _stressType, PL_NUM _p0, PL_NUM _tauP );
	void setSwitchTime( N_PRES _switchTime );
	void setResArray( N_PRES* _resArr );
	void setResDtArray( N_PRES* _resArrDt );

	PL_NUM increaseTime();

	void calc_nonlin_system_run_test( long  _x, long _t );

	PL_NUM do_step();
	void dump_sol();
	void dump_check_sol( int fNum );
	void dumpSolAll( int fNum );
	void dump_left_border_vals();
	void dumpMatrA();
	void dump_whole_sol( int var );

	N_PRES cur_t;
	N_PRES dt;			//time step
	N_PRES switchTime;
	int curTimeStep;

	SolverPar saveParamsToStruct();

	time_t getOrthoBTime();
	time_t totalTime;
	time_t totalTime1;
	time_t rgkTime;
	time_t rgkTime1;
	time_t orthoTime;
	time_t orthoTime1;
	time_t matrTime;
	time_t matrTime1;
	time_t buildSolnTime;
	time_t buildSolnTime1;
private:
	void calcConsts();

	int maxNewtonIterReached;

	PL_NUM E1;				//Young's modulus
	PL_NUM E2;				//Young's modulus
	PL_NUM nu21;			//Poisson's ratio	
	PL_NUM rho;				//composite's density

	PL_NUM sigma_x;			//electric conductivity
	PL_NUM sigma_x_mu;

	N_PRES h;				//thickness of the plate
	N_PRES a;				//width of the plate

	int eq_num;				//number of equations in main system

	PL_NUM J0;
	PL_NUM J0_1;
	PL_NUM  p0;				//constant mechanical load
	int stressType;
	int currentType;

	PL_NUM B11;
	PL_NUM B22;
	PL_NUM B12;
	PL_NUM By0;
	PL_NUM By1;                                      // in considered boundary-value problem
	PL_NUM By2;										//CAUTION! almost everywhere By2 is assumed to be 0. If that is not true, the program must be fixed

	PL_NUM eps_0;
	PL_NUM eps_x;
	PL_NUM eps_x_0;

	PL_NUM tauP;
	//PL_NUM tauJ;
	PL_NUM tauSin;
	PL_NUM tauSin_1;
	PL_NUM tauExp;
	PL_NUM tauExp_1;
	PL_NUM omega;
	PL_NUM omega_1;
	PL_NUM rad;

	int Km;				//number of steps by space
	//int Kt;				//number of steps by time

	N_PRES dx;

	N_PRES al;			//some weird var for normalization. It is said that this will improve the numerical scheme. must be equal to density
	PL_NUM beta;		//parameter at Newmark's time integration scheme

	vector<VarVect<PL_NUM> > mesh;		//2d mesh for solution.
	//vector<PL_NUM> nonlin_matr_A;		//matrix A for the nonlinear system at certain t and x
	//vector<PL_NUM> nonlin_vect_f;		//vector f on right part of nonlinear system at certain t and x
	PL_NUM nonlin_matr_A[EQ_NUM][EQ_NUM];		//matrix A for the nonlinear system at certain t and x

	//PL_NUM nonlin_vect_f[EQ_NUM];		//vector f on right part of nonlinear system at certain t and x

	//Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor> nonlin_matr_A;
	Matrix<PL_NUM, EQ_NUM, 1> nonlin_vect_f;
	
	Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor> matrA;
	Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor> matrA1;
	Matrix<PL_NUM, EQ_NUM, 1> vectF;
	Matrix<PL_NUM, EQ_NUM, 1> vectF1;

	vector<PL_NUM> newmark_A;
	vector<PL_NUM> newmark_B;

	Matrix<PL_NUM, EQ_NUM, 1> N1;
	Matrix<PL_NUM, EQ_NUM, 1> N2;
	Matrix<PL_NUM, EQ_NUM, 1> N3;
	Matrix<PL_NUM, EQ_NUM, 1> N4;
	Matrix<PL_NUM, EQ_NUM, 1> N5;

	Matrix<PL_NUM, EQ_NUM, 1> N1orthog;
	Matrix<PL_NUM, EQ_NUM, 1> N2orthog;
	Matrix<PL_NUM, EQ_NUM, 1> N3orthog;
	Matrix<PL_NUM, EQ_NUM, 1> N4orthog;
	Matrix<PL_NUM, EQ_NUM, 1> N5orthog;

	Matrix<PL_NUM, EQ_NUM, 1> tmpPhi1;
	Matrix<PL_NUM, EQ_NUM, 1> tmpPhi2;
	Matrix<PL_NUM, EQ_NUM, 1> tmpPhi3;
	Matrix<PL_NUM, EQ_NUM, 1> tmpPhi4;
	Matrix<PL_NUM, EQ_NUM, 1> tmpPhi5;

	PL_NUM Phi1[ABM_STAGE_NUM][EQ_NUM];
	PL_NUM Phi2[ABM_STAGE_NUM][EQ_NUM];
	PL_NUM Phi3[ABM_STAGE_NUM][EQ_NUM];
	PL_NUM Phi4[ABM_STAGE_NUM][EQ_NUM];
	PL_NUM Phi5[ABM_STAGE_NUM][EQ_NUM];

	RungeKutta<PL_NUM>* rungeKutta;
	OrthoBuilder<PL_NUM>* orthoBuilder;

	void calc_Newmark_AB( int _x );		//don't really know why we need this stuff with mode need to FIX -- UPD "mode" removed
	void calc_nonlin_system( int _x );	//though in fact, the system is linear, "nonlinear" means that this is for the solution of the nonlinear problem
	void calcLinSystem( int _x, Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor>* A, Matrix<PL_NUM, EQ_NUM, 1>* f );

	int checkConv();
	void copyToResArr();
};

template<class PL_NUM>
Solver<PL_NUM>::Solver() :
	rungeKutta( 0 ),
	orthoBuilder( 0 ),
	totalTime1( 0 ),
	totalTime( 0 ),
	rgkTime1( 0 ),
	rgkTime( 0 ),
	orthoTime1( 0 ),
	orthoTime( 0 ),
	matrTime1( 0 ),
	matrTime( 0 ),
	buildSolnTime1( 0 ),
	buildSolnTime( 0 ),
	resArr( 0 ),
	resArrDt( 0 ),
	maxNewtonIterReached( 0 )
{
	
}

template<class PL_NUM>
int Solver<PL_NUM>::getMaxNewtonIterReached()
{
	return maxNewtonIterReached;
}

template<class PL_NUM>
void Solver<PL_NUM>::copyToResArr() 	//do nothing in this case
{
	//cout << " -- Solver does not keep the solution\n";
}

template<>
void inline Solver<N_PRES>::copyToResArr() 	//do nothing in this case
{
	if( resArr != 0 && resArrDt != 0 )
	{
		for( int y = 0; y < Km; ++y )
		{
			for( int i = 0; i < eq_num; ++i )
			{
				resArr[curTimeStep * ( Km * eq_num ) + y * eq_num + i] = mesh[y].Nk1[i];
				resArrDt[curTimeStep * ( Km * eq_num ) + y * eq_num + i] = mesh[y].d1N0[i];
			}
		}
	}
}

template<class PL_NUM>
SolverPar Solver<PL_NUM>::saveParamsToStruct() 	//do nothing in this case
{
	cout << " -- Solver does not save the params\n";
}

template<>
SolverPar inline Solver<N_PRES>::saveParamsToStruct()
{
	SolverPar saveTo;

	saveTo.E1 = E1;
	saveTo.E2 = E2;
	saveTo.nu21 = nu21;
	saveTo.rho = rho;
	
	saveTo.sigma_x = sigma_x;
	saveTo.sigma_x_mu = sigma_x_mu;

	saveTo.h = h;
	saveTo.a = a;

	saveTo.dt = dt;
	saveTo.Km = Km;
	saveTo.dx = dx;

	saveTo.J0 = J0;
	saveTo.J0_1 = J0_1;
	saveTo.tauSin = tauSin;
	saveTo.tauSin_1 = tauSin_1;
	saveTo.tauExp = tauExp;
	saveTo.tauExp_1 = tauExp_1;

	saveTo.p0 = p0;				//constant mechanical load
	saveTo.tauP = tauP;
	saveTo.rad = rad;

	saveTo.stressType = stressType;
	saveTo.currentType = currentType;

	saveTo.beta = beta;

	saveTo.B11 = B11;
	saveTo.B22 = B22;
	saveTo.B12 = B12;
	saveTo.By0 = By0;
	saveTo.By1 = By1;                                      // in considered boundary-value problem

	saveTo.eps_0 = eps_0;
	saveTo.eps_x = eps_x;
	saveTo.eps_x_0 = eps_x_0;

	return saveTo;
}

template<class PL_NUM>
Solver<PL_NUM>::~Solver()
{
	delete rungeKutta;
	delete orthoBuilder;
}

template<class PL_NUM>
void Solver<PL_NUM>::setResArray( N_PRES* _resArr )
{
	if( _resArr != 0 )
	{
		resArr = _resArr;
		for( int i = 0; i < ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM; ++i )
		{
			resArr[i] = 0.0;
		}
	}
	else
	{
		cout << "WARNING! pointer to the solution of the primal problem is NULL!\n";
	}
}

template<class PL_NUM>
void Solver<PL_NUM>::setResDtArray( N_PRES* _resArrDt )
{
	if( _resArrDt != 0 )
	{
		resArrDt = _resArrDt;
		for( int i = 0; i < ( int )( CHAR_TIME / DELTA_T  + 1 ) * NODES_Y * EQ_NUM; ++i )
		{
			resArrDt[i] = 0.0;
		}
	}
	else
	{
		cout << "WARNING! pointer to the solution of the dt of the primal problem is NULL!\n";
	}
}

template<class PL_NUM>
time_t Solver<PL_NUM>::getOrthoBTime()
{
	return orthoBuilder->orthoTotal;
}

template<class PL_NUM>
void Solver<PL_NUM>::setMechLoad( int _stressType, PL_NUM _p0, PL_NUM _tauP )
{
	stressType = _stressType;
	p0 = _p0;
	tauP = _tauP;
	rad = 0.0021 / 100.0;
}

template<class PL_NUM>
void Solver<PL_NUM>::setSwitchTime( N_PRES _switchTime )
{
	switchTime = _switchTime;
}

template<class PL_NUM>
void Solver<PL_NUM>::setTask( PL_NUM _J0, PL_NUM _tauSin, PL_NUM _tauExp,
							PL_NUM _J0_1, PL_NUM _tauSin_1, PL_NUM _tauExp_1,
							PL_NUM _By0, int _stressType, PL_NUM _p0, PL_NUM _tauP )
{
	totalTime1 = 0;
	totalTime = 0;
	rgkTime1 = 0;
	rgkTime = 0;
	orthoTime1 = 0;
	orthoTime = 0;

	ofstream of1( "test_sol.txt" );
	of1.close();

//material properties
	al = 1.0;
	E1 = 102970000000.0;			//Young's modulus
	E2 = 7550000000.0;				//Young's modulus
	nu21 = 0.3;						//Poisson's ratio	
	rho = 1594.0;					//composite's density

	sigma_x = 39000.0;				//electric conductivity
	sigma_x_mu = sigma_x * 0.000001256l;

	h = 0.0021;						//thickness of the plate
	a = 0.1524;						//width of the plate
//other
	tauP = _tauP;//0.01;
	p0 = _p0;//10000000;
	rad = h / 100.0;

	_tauSin != 0.0l ? tauSin = _tauSin : tauSin = 1;
	_tauExp != 0.0l ? tauExp = _tauExp : tauExp = 1;
	_tauSin_1 != 0.0l ? tauSin_1 = _tauSin_1 : tauSin_1 = 1;
	_tauExp_1 != 0.0l ? tauExp_1 = _tauExp_1 : tauExp_1 = 1;

	eq_num = EQ_NUM;

	J0 = _J0;
	J0 *= J0_SCALE;
	J0_1 = _J0_1;
	J0_1 *= J0_SCALE;

	omega = (long double)M_PI / tauSin;
	omega_1 = (long double)M_PI / tauSin_1;
	
	stressType = _stressType;
	currentType = current_exp_sin;

	By0 = _By0;
	By0 *= BY0_SCALE;

	eps_0 = 0.000000000008854;
	eps_x = 0.0000000002501502912;

	Km = NODES_Y;
	dx = a / ( Km - 1 );

	cout << " dx and rad are: " << dx << " " << rad << endl;

	dt = DELTA_T;
	cur_t = dt;
	curTimeStep = 1;
	switchTime = SWITCH_TIME;

	beta = 0.25;

	if( rungeKutta != 0 )
	{
		delete rungeKutta;
	}
	if( orthoBuilder != 0 )
	{
		delete orthoBuilder;
	}

	rungeKutta = new RungeKutta<PL_NUM>( eq_num );
	orthoBuilder = new OrthoBuilderGSh<PL_NUM>();
	orthoBuilder->setParams( Km );

	mesh.resize( 0 );
	mesh.resize( Km );
	for( int i = 0; i < mesh.size(); ++i ){
		mesh[i].setup( eq_num );
	}
	//solInfoMap.resize( Km );

	for( int i = 0; i < EQ_NUM; ++i )
	{
		for( int j = 0; j < EQ_NUM; ++j )
		{
			nonlin_matr_A[i][j] = 0.0l;
			matrA( i, j ) = 0.0l;
			matrA1( i, j ) = 0.0l;
		}
		nonlin_vect_f( i ) = 0.0l;
		vectF( i ) = 0.0l;
		vectF1( i ) = 0.0l;
	}

	newmark_A.resize( eq_num );
	newmark_B.resize( eq_num );

	for( int i = 0; i < EQ_NUM; ++i )
	{
		newmark_A[i] = 0.0l;
		newmark_B[i] = 0.0l;
		N1( i ) = 0.0l;
		N2( i ) = 0.0l;
		N3( i ) = 0.0l;
		N4( i ) = 0.0l;
		N5( i ) = 0.0l;
	}

	for( int i = 0; i < ABM_STAGE_NUM; ++i )
	{
		for( int j = 0; j < EQ_NUM; ++j )
		{
			Phi1[i][j] = 0.0l;
			Phi2[i][j] = 0.0l;
			Phi3[i][j] = 0.0l;
			Phi4[i][j] = 0.0l;
			Phi5[i][j] = 0.0l;
		}
	}

	calcConsts();
}

template<class PL_NUM>
void Solver<PL_NUM>::calcConsts()
{
	B11 = E1 * E1 / ( E1 - nu21 * nu21 * E2 );
	B22 = E2 / ( 1.0l - nu21 * nu21 * E2 / E1 );
	B12 = nu21 * E2 * E1 / ( E1 - nu21 * nu21 * E2 );

	By1 = 2.0l * By0;                                      // in considered boundary-value problem
	By2 = 0.0l;
	eps_x_0 = eps_x - eps_0;
}

template<class PL_NUM>
PL_NUM Solver<PL_NUM>::increaseTime()
{
	cur_t += dt;
	++curTimeStep;
	return cur_t;
}

template<class PL_NUM>
void Solver<PL_NUM>::calc_Newmark_AB( int _x )
{
	for( int i = 0; i < eq_num; ++i)
	{
		newmark_A[i] = -mesh[_x].Nk0[i] / beta / dt / dt - mesh[_x].d1N0[i] / beta / dt
						- ( 0.5l - beta ) / beta * mesh[_x].d2N0[i];
		newmark_B[i] = -0.5l * mesh[_x].Nk0[i] / beta / dt + ( 1.0l - 0.5l / beta ) * mesh[_x].d1N0[i]
						- 0.5l * dt * ( 1.0l- 1.0l / beta * ( 0.5l - beta ) ) * mesh[_x].d2N0[i];
	}
}

template<class PL_NUM>
void Solver<PL_NUM>::calc_nonlin_system( int _x )
{
	PL_NUM Jx = 0.0;
	if( currentType == current_const )
	{
		Jx = J0;
	}
	else if( currentType == current_sin )
	{
		Jx = J0 * sin( omega * cur_t );
	}
	else if( currentType == current_exp_sin )
	{
		if( cur_t <= switchTime )
		{
			Jx = J0 * exp( -cur_t / tauExp ) * sin( omega * cur_t );
		}
		else
		{
			Jx = J0_1 * exp( -cur_t / tauExp_1 ) * sin( omega_1 * cur_t );
		}
	}
	PL_NUM Pimp = 0.0l;
	if( stressType == stress_centered )
	{
		if( cur_t < tauP && fabs( (long double)_x * dx - a / 2.0 ) < rad )
		{
			Pimp = p0 * sqrt( 1.0l - fabs( (long double)_x * dx - a / 2.0l ) * fabs( (long double)_x * dx - a / 2.0 ) / rad / rad	) 
				* sin( (long double)M_PI * cur_t / tauP );
		}
	}
	else if( stressType == stress_whole )
	{
		Pimp = p0;
	}

	nonlin_matr_A[ 0][3 ] = 1.0l / al / h / B22;
	nonlin_matr_A[ 1][2 ] = 1.0l / al;
	nonlin_matr_A[ 2][5 ] = -12.0l /al / h / h / h / B22;

	nonlin_matr_A[ 3][0 ] = rho / al * h / beta / dt / dt + sigma_x * h / 2.0l / beta / dt / al * mesh[_x].Nk[7] * mesh[_x].Nk[7];
	nonlin_matr_A[ 3][1 ] = -sigma_x * h / 4.0l / beta / dt / al * By1 * mesh[_x].Nk[7];
	nonlin_matr_A[ 3][2 ] = -eps_x_0 / 4.0l / beta / dt / al * h * By1 * mesh[_x].Nk[6];
	nonlin_matr_A[ 3][3 ] = eps_x_0 / 2.0l / beta / dt / B22 / al * h * mesh[_x].Nk[6] * mesh[_x].Nk[7];
	nonlin_matr_A[ 3][6 ] = 1.0l / al * ( sigma_x * h * mesh[_x].Nk[7] + eps_x_0 * h / B22 * mesh[_x].Nk[7] * ( 1.0 / 2.0l / beta / dt 
									* mesh[_x].Nk[3]
									+ newmark_B[3] ) - eps_x_0 / 2.0l * h * By1 * ( 1.0l / 2.0l / beta / dt * mesh[_x].Nk[2] + newmark_B[2] ) );
	nonlin_matr_A[ 3][7 ] = 1.0l / al * ( sigma_x * h * mesh[_x].Nk[6] + 2.0l * sigma_x * h * mesh[_x].Nk[7] 
									* ( 1.0 / 2.0l / beta / dt * mesh[_x].Nk[0] + newmark_B[0] )
									- sigma_x * h / 2.0l * By1 * ( 1.0l / 2.0l / beta / dt * mesh[_x].Nk[1] + newmark_B[1] ) 
									+ eps_x_0 * h / B22 * mesh[_x].Nk[6]
									* ( 1.0 / 2.0l / beta / dt * mesh[_x].Nk[3] 
									+ newmark_B[3] ) + h * Jx );

	nonlin_matr_A[ 4][0 ] = -sigma_x * h / 4.0l / beta / dt / al * By1 * mesh[_x].Nk[7];
	nonlin_matr_A[ 4][1 ] = rho / al * h / beta / dt / dt + sigma_x * h / 8.0l / beta / al / dt * ( By1 * By1 );
	nonlin_matr_A[ 4][2 ] = 1.0 / 2.0l / beta / al * ( -eps_x_0 * h * mesh[_x].Nk[6] * mesh[_x].Nk[7] ) / dt;
	nonlin_matr_A[ 4][6 ] = -1.0l / al * ( sigma_x * h / 2.0l * By1 + eps_x_0 * h 
									* mesh[_x].Nk[7] * ( 1.0 / 2.0l / beta * mesh[_x].Nk[2] / dt + newmark_B[2] ) );
	nonlin_matr_A[ 4][7 ] = -1.0l / al * ( sigma_x * h / 2.0l * By1 * ( 1.0 / 2.0l / beta * mesh[_x].Nk[0] / dt 
									+ newmark_B[0] ) - ( - eps_x_0 * h * mesh[_x].Nk[6] ) * ( 1.0 / 2.0l / beta * mesh[_x].Nk[2] / dt + newmark_B[2] ) );

	//nonlin_matr_A[ 5][1 ] = -sigma_x * h * h / 24.0l / beta / dt / al * By2 * mesh[_x].Nk[7];	//because By2 is always 0 for the problems that we consider
	nonlin_matr_A[ 5][2 ] = -1.0 / 2.0l / beta / dt / al * ( sigma_x * h * h * h / 12.0l * mesh[_x].Nk[7] 
									* mesh[_x].Nk[7] ) - h * h * h / 12.0l / beta / dt / dt * rho / al;
	nonlin_matr_A[ 5][4 ] = 1.0l / al;
	nonlin_matr_A[ 5][5 ] = eps_x_0 / 2.0l / beta / B22 / dt / al * mesh[_x].Nk[6] * mesh[_x].Nk[7];
	nonlin_matr_A[ 5][6 ] = -eps_x_0 / al * ( -mesh[_x].Nk[7] / B22 * ( 1.0l / 2.0l / beta / dt * mesh[_x].Nk[5] + newmark_B[5] ) );
	nonlin_matr_A[ 5][7 ] = -1.0l / al * ( sigma_x * h * h * h / 6.0l * mesh[_x].Nk[7] * ( 1.0l / 2.0l / beta / dt * mesh[_x].Nk[2] 
									+ newmark_B[2] ) - eps_x_0 / B22 * mesh[_x].Nk[6] * ( 1.0l / 2.0l / beta / dt * mesh[_x].Nk[5] + newmark_B[5] ) );

	nonlin_matr_A[ 6][7 ] = 1.0l / 2.0l / ( beta * al * dt );

	nonlin_matr_A[ 7][0 ] = sigma_x_mu / 2.0l / beta / ( dt * al ) * mesh[_x].Nk[7];
	nonlin_matr_A[ 7][1 ] = -sigma_x_mu / 4.0l / beta / ( dt * al ) * By1;
	nonlin_matr_A[ 7][6 ] = sigma_x_mu / al;
	nonlin_matr_A[ 7][7 ] = sigma_x_mu / al * ( 1.0 / 2.0l / beta / dt * mesh[_x].Nk[0] + newmark_B[0]);

	nonlin_vect_f( 3 ) = rho / al * h * newmark_A[0] + 1.0l / al * ( -sigma_x * h * mesh[_x].Nk[6] 
						* mesh[_x].Nk[7] - sigma_x * h / beta / dt * mesh[_x].Nk[7] 
						* mesh[_x].Nk[7] * mesh[_x].Nk[0] - sigma_x * h * mesh[_x].Nk[7] 
						* mesh[_x].Nk[7] * newmark_B[0] 
						+ sigma_x * h / 4.0l / beta / dt * By1 * mesh[_x].Nk[7] * mesh[_x].Nk[1] 
						- eps_x_0 * h / beta / dt / B22 * mesh[_x].Nk[6] * mesh[_x].Nk[7] * mesh[_x].Nk[3] - eps_x_0 * h / B22 
						* mesh[_x].Nk[6] * mesh[_x].Nk[7] * newmark_B[3] 
						+ eps_x_0 * h / 4.0l / beta / dt * By1 * mesh[_x].Nk[6] * mesh[_x].Nk[2] );
	nonlin_vect_f( 4 ) = rho / al * h * newmark_A[1] + Pimp / al + 1.0l / al * ( sigma_x * h / 4.0l / beta * By1 / dt	//55454 h
						* mesh[_x].Nk[7] * mesh[_x].Nk[0] + sigma_x * h / 4.0l * ( By1 * By1 ) * newmark_B[1] 
						+ eps_x_0 * h / beta / dt 
						* mesh[_x].Nk[6] * mesh[_x].Nk[7] * mesh[_x].Nk[2] + eps_x_0 * h * mesh[_x].Nk[6] 
						* mesh[_x].Nk[7] * newmark_B[2] - h / 2.0l * Jx * By1 );			//By1
	nonlin_vect_f( 5 ) = 1.0l / al * ( sigma_x * h * h * h / 12.0l / beta / dt * mesh[_x].Nk[7] * mesh[_x].Nk[7] * mesh[_x].Nk[2] 
						+ sigma_x * h * h * h / 12.0l * mesh[_x].Nk[7] * mesh[_x].Nk[7] * newmark_B[2] - eps_x_0 / beta / dt / B22 * mesh[_x].Nk[6] 
						* mesh[_x].Nk[7] * mesh[_x].Nk[5] - eps_x_0 / B22 * mesh[_x].Nk[6] * mesh[_x].Nk[7] 
						* newmark_B[5] ) - h * h * h / 12.0l * newmark_A[2] * rho  / al;
	nonlin_vect_f( 6 ) = newmark_B[7] / al;
	nonlin_vect_f( 7 ) = 1.0l / al * ( -sigma_x_mu / 2.0l / beta / dt * mesh[_x].Nk[7] * mesh[_x].Nk[0] - 0.5l * sigma_x_mu * By1 * newmark_B[1] );
}

template<class PL_NUM>
void Solver<PL_NUM>::calcLinSystem( int _x, Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor>* A, Matrix<PL_NUM, EQ_NUM, 1>* f )
{
	PL_NUM Jx = 0.0;
	if( currentType == current_const )
	{
		Jx = J0;
	}
	else if( currentType == current_sin )
	{
		Jx = J0 * sin( omega * cur_t );
	}
	else if( currentType == current_exp_sin )
	{
		if( cur_t <= switchTime )
		{
			Jx = J0 * exp( -cur_t / tauExp ) * sin( omega * cur_t );
		}
		else
		{
			Jx = J0_1 * exp( -cur_t / tauExp_1 ) * sin( omega_1 * cur_t );
		}
	}
	PL_NUM Pimp = 0.0l;
	if( stressType == stress_centered )
	{
		if( cur_t < tauP && fabs( (long double)_x * dx - a / 2.0 ) < rad )
		{
			Pimp = p0 * sqrt( 1.0l - fabs( (long double)_x * dx - a / 2.0l ) * fabs( (long double)_x * dx - a / 2.0 ) / rad / rad	) 
				* sin( (long double)M_PI * cur_t / tauP );
		}
	}
	else if( stressType == stress_whole )
	{
		Pimp = p0;
	}

	(*A)( 0, 3 ) = 1.0l / al / h / B22;
	(*A)( 1, 2 ) = 1.0l / al;
	(*A)( 2, 5 ) = -12.0l /al / h / h / h / B22;

	(*A)( 3, 0 ) = rho / al * h / beta / dt / dt + sigma_x * h / 2.0l / beta / dt / al * mesh[_x].Nk[7] * mesh[_x].Nk[7];
	(*A)( 3, 1 ) = -sigma_x * h / 4.0l / beta / dt / al * By1 * mesh[_x].Nk[7];
	(*A)( 3, 2 ) = -eps_x_0 / 4.0l / beta / dt / al * h * By1 * mesh[_x].Nk[6];
	(*A)( 3, 3 ) = eps_x_0 / 2.0l / beta / dt / B22 / al * h * mesh[_x].Nk[6] * mesh[_x].Nk[7];
	(*A)( 3, 6 ) = 1.0l / al * ( sigma_x * h * mesh[_x].Nk[7] + eps_x_0 * h / B22 * mesh[_x].Nk[7] * ( 1.0 / 2.0l / beta / dt 
									* mesh[_x].Nk[3]
									+ newmark_B[3] ) - eps_x_0 / 2.0l * h * By1 * ( 1.0l / 2.0l / beta / dt * mesh[_x].Nk[2] + newmark_B[2] ) );
	(*A)( 3, 7 ) = 1.0l / al * ( sigma_x * h * mesh[_x].Nk[6] + 2.0l * sigma_x * h * mesh[_x].Nk[7] 
									* ( 1.0 / 2.0l / beta / dt * mesh[_x].Nk[0] + newmark_B[0] )
									- sigma_x * h / 2.0l * By1 * ( 1.0l / 2.0l / beta / dt * mesh[_x].Nk[1] + newmark_B[1] ) 
									+ eps_x_0 * h / B22 * mesh[_x].Nk[6]
									* ( 1.0 / 2.0l / beta / dt * mesh[_x].Nk[3] 
									+ newmark_B[3] ) + h * Jx );

	(*A)( 4, 0 ) = -sigma_x * h / 4.0l / beta / dt / al * By1 * mesh[_x].Nk[7];
	(*A)( 4, 1 ) = rho / al * h / beta / dt / dt + sigma_x * h / 8.0l / beta / al / dt * ( By1 * By1 );
	(*A)( 4, 2 ) = 1.0 / 2.0l / beta / al * ( -eps_x_0 * h * mesh[_x].Nk[6] * mesh[_x].Nk[7] ) / dt;
	(*A)( 4, 6 ) = -1.0l / al * ( sigma_x * h / 2.0l * By1 + eps_x_0 * h 
									* mesh[_x].Nk[7] * ( 1.0 / 2.0l / beta * mesh[_x].Nk[2] / dt + newmark_B[2] ) );
	(*A)( 4, 7 ) = -1.0l / al * ( sigma_x * h / 2.0l * By1 * ( 1.0 / 2.0l / beta * mesh[_x].Nk[0] / dt 
									+ newmark_B[0] ) - ( - eps_x_0 * h * mesh[_x].Nk[6] ) * ( 1.0 / 2.0l / beta * mesh[_x].Nk[2] / dt + newmark_B[2] ) );

	(*A)( 5, 2 ) = -1.0 / 2.0l / beta / dt / al * ( sigma_x * h * h * h / 12.0l * mesh[_x].Nk[7] 
									* mesh[_x].Nk[7] ) - h * h * h / 12.0l / beta / dt / dt * rho / al;
	(*A)( 5, 4 ) = 1.0l / al;
	(*A)( 5, 5 ) = eps_x_0 / 2.0l / beta / B22 / dt / al * mesh[_x].Nk[6] * mesh[_x].Nk[7];
	(*A)( 5, 6 ) = -eps_x_0 / al * ( -mesh[_x].Nk[7] / B22 * ( 1.0l / 2.0l / beta / dt * mesh[_x].Nk[5] + newmark_B[5] ) );
	(*A)( 5, 7 ) = -1.0l / al * ( sigma_x * h * h * h / 6.0l * mesh[_x].Nk[7] * ( 1.0l / 2.0l / beta / dt * mesh[_x].Nk[2] 
									+ newmark_B[2] ) - eps_x_0 / B22 * mesh[_x].Nk[6] * ( 1.0l / 2.0l / beta / dt * mesh[_x].Nk[5] + newmark_B[5] ) );

	(*A)( 6, 7 ) = 1.0l / 2.0l / ( beta * al * dt );

	(*A)( 7, 0 ) = sigma_x_mu / 2.0l / beta / ( dt * al ) * mesh[_x].Nk[7];
	(*A)( 7, 1 ) = -sigma_x_mu / 4.0l / beta / ( dt * al ) * By1;
	(*A)( 7, 6 ) = sigma_x_mu / al;
	(*A)( 7, 7 ) = sigma_x_mu / al * ( 1.0 / 2.0l / beta / dt * mesh[_x].Nk[0] + newmark_B[0]);

	(*f)( 3 ) = rho / al * h * newmark_A[0] + 1.0l / al * ( -sigma_x * h * mesh[_x].Nk[6] 
						* mesh[_x].Nk[7] - sigma_x * h / beta / dt * mesh[_x].Nk[7] 
						* mesh[_x].Nk[7] * mesh[_x].Nk[0] - sigma_x * h * mesh[_x].Nk[7] 
						* mesh[_x].Nk[7] * newmark_B[0] 
						+ sigma_x * h / 4.0l / beta / dt * By1 * mesh[_x].Nk[7] * mesh[_x].Nk[1] 
						- eps_x_0 * h / beta / dt / B22 * mesh[_x].Nk[6] * mesh[_x].Nk[7] * mesh[_x].Nk[3] - eps_x_0 * h / B22 
						* mesh[_x].Nk[6] * mesh[_x].Nk[7] * newmark_B[3] 
						+ eps_x_0 * h / 4.0l / beta / dt * By1 * mesh[_x].Nk[6] * mesh[_x].Nk[2] );
	(*f)( 4 ) = rho / al * h * newmark_A[1] + Pimp / al + 1.0l / al * ( sigma_x * h / 4.0l / beta * By1 / dt	//55454 h
						* mesh[_x].Nk[7] * mesh[_x].Nk[0] + sigma_x * h / 4.0l * ( By1 * By1 ) * newmark_B[1] 
						+ eps_x_0 * h / beta / dt 
						* mesh[_x].Nk[6] * mesh[_x].Nk[7] * mesh[_x].Nk[2] + eps_x_0 * h * mesh[_x].Nk[6] 
						* mesh[_x].Nk[7] * newmark_B[2] - h / 2.0l * Jx * By1 );			//By1
	(*f)( 5 ) = 1.0l / al * ( sigma_x * h * h * h / 12.0l / beta / dt * mesh[_x].Nk[7] * mesh[_x].Nk[7] * mesh[_x].Nk[2] 
						+ sigma_x * h * h * h / 12.0l * mesh[_x].Nk[7] * mesh[_x].Nk[7] * newmark_B[2] - eps_x_0 / beta / dt / B22 * mesh[_x].Nk[6] 
						* mesh[_x].Nk[7] * mesh[_x].Nk[5] - eps_x_0 / B22 * mesh[_x].Nk[6] * mesh[_x].Nk[7] 
						* newmark_B[5] ) - h * h * h / 12.0l * newmark_A[2] * rho  / al;
	(*f)( 6 ) = newmark_B[7] / al;
	(*f)( 7 ) = 1.0l / al * ( -sigma_x_mu / 2.0l / beta / dt * mesh[_x].Nk[7] * mesh[_x].Nk[0] - 0.5l * sigma_x_mu * By1 * newmark_B[1] );
}

template<class PL_NUM>
PL_NUM Solver<PL_NUM>::do_step()
{	
	totalTime1 = time( 0 );
	//cout << "cur time is " << cur_t.real() << endl;
	//cout << "time step number " << curTimeStep << endl;
	//cout << "calculating solution for the next time step\n\n";

	int iter( 0 );
	int preLin( 0 );
	int cont( 1 );

	do{
		//cout << " proc num " << omp_get_thread_num() << endl;

		int active = 1;

		for( int x = 0; x < Km; ++x )
		{
			for( int i = 0; i < eq_num; ++i )
			{
				mesh[x].Nk[i] = mesh[x].Nk1[i];
			}
		}

		if( iter == 0 && curTimeStep == 1 )
		{
			preLin = 0;
		}
		else
		{
			preLin = 1;
		}
		calc_Newmark_AB( 0 );
		calcLinSystem( 0, &matrA1, &vectF1 );

		orthoBuilder->flushO( 0 );
		orthoBuilder->resetOrthoDoneInfo();

		N1( 0 ) = 0.0; N2( 0 ) = 0.0; N3( 0 ) = 0.0; N4( 0 ) = 0.0;     
		N1( 1 ) = 0.0; N2( 1 ) = 0.0; N3( 1 ) = 0.0; N4( 1 ) = 0.0;
		N1( 2 ) = 1.0; N2( 2 ) = 0.0; N3( 2 ) = 0.0; N4( 2 ) = 0.0;
		N1( 3 ) = 0.0; N2( 3 ) = 1.0; N3( 3 ) = 0.0; N4( 3 ) = 0.0;
		N1( 4 ) = 0.0; N2( 4 ) = 0.0; N3( 4 ) = 1.0; N4( 4 ) = 0.0;
		N1( 5 ) = 0.0; N2( 5 ) = 0.0; N3( 5 ) = 0.0; N4( 5 ) = 0.0; 
		N1( 6 ) = 0.0; N2( 6 ) = 0.0; N3( 6 ) = 0.0; N4( 6 ) = mesh[0].Nk[0] / 2.0l / beta / dt + newmark_B[0];
		N1( 7 ) = 0.0; N2( 7 ) = 0.0; N3( 7 ) = 0.0; N4( 7 ) = -1.0;

		N5( 0 ) = 0.0; 
		N5( 1 ) = 0.0; 
		N5( 2 ) = 0.0; 
		N5( 3 ) = 0.0; 
		N5( 4 ) = 0.0; 
		N5( 5 ) = 0.0; 
		N5( 6 ) = -( newmark_B[0] * mesh[0].Nk[7] - newmark_B[1] * By0 );
		N5( 7 ) = mesh[0].Nk[7];

		orthoBuilder->setInitVects( N1, N2, N3, N4, N5 );

		for( int x = 0; x < Km - 1; ++x )
		{
			orthoBuilder->flushO( x + 1 );

			//calc_Newmark_AB( x );

			//matrTime1 = time( 0 );
			//calc_nonlin_system( x );
			//matrTime += time( 0 ) - matrTime1;

			matrA = matrA1;
			vectF = vectF1;

			tmpPhi1 = matrA * N1;
			tmpPhi2 = matrA * N2;
			tmpPhi3 = matrA * N3;
			tmpPhi4 = matrA * N4;
			tmpPhi5 = matrA * N5 + vectF;

			int PhiInd = 0;

			if( active <= ABM_STAGE_NUM )
			{
				PhiInd = active - 1;
			}
			else
			{
				for( int i = 0; i < ABM_STAGE_NUM - 1; ++i )
				{
					for( int j = 0; j < EQ_NUM; ++j )
					{
						Phi1[i][j] = Phi1[i + 1][j];
						Phi2[i][j] = Phi2[i + 1][j];
						Phi3[i][j] = Phi3[i + 1][j];
						Phi4[i][j] = Phi4[i + 1][j];
						Phi5[i][j] = Phi5[i + 1][j];
					}
				}
				PhiInd = ABM_STAGE_NUM - 1;
			}
		
			for( int i = 0; i < EQ_NUM; ++i )
			{
				Phi1[PhiInd][i] = tmpPhi1( i );
				Phi2[PhiInd][i] = tmpPhi2( i );
				Phi3[PhiInd][i] = tmpPhi3( i );
				Phi4[PhiInd][i] = tmpPhi4( i );
				Phi5[PhiInd][i] = tmpPhi5( i );
			}

			calc_Newmark_AB( x + 1 );
			calcLinSystem( x + 1, &matrA1, &vectF1 );

			rgkTime1 = time( 0 );

			/*rungeKutta->calc( nonlin_matr_A, nonlin_vect_f, dx, 0, &N1 );
			rungeKutta->calc( nonlin_matr_A, nonlin_vect_f, dx, 0, &N2 );
			rungeKutta->calc( nonlin_matr_A, nonlin_vect_f, dx, 0, &N3 );
			rungeKutta->calc( nonlin_matr_A, nonlin_vect_f, dx, 0, &N4 );
			rungeKutta->calc( nonlin_matr_A, nonlin_vect_f, dx, 1, &N5 );*/

			if( active >= 4 )
			{
				//use ABM method
				//predictor:
				for( int i = 0; i < EQ_NUM; ++i )
				{
					tmpPhi1( i ) = N1( i ) + dx / 24.0l
									* ( 55.0l * Phi1[3][i] - 59.0l * Phi1[2][i] + 37.0l * Phi1[1][i] - 9.0l * Phi1[0][i] );
					tmpPhi2( i ) = N2( i ) + dx / 24.0l
									* ( 55.0l * Phi2[3][i] - 59.0l * Phi2[2][i] + 37.0l * Phi2[1][i] - 9.0l * Phi2[0][i] );
					tmpPhi3( i ) = N3( i ) + dx / 24.0l
									* ( 55.0l * Phi3[3][i] - 59.0l * Phi3[2][i] + 37.0l * Phi3[1][i] - 9.0l * Phi3[0][i] );
					tmpPhi4( i ) = N4( i ) + dx / 24.0l
									* ( 55.0l * Phi4[3][i] - 59.0l * Phi4[2][i] + 37.0l * Phi4[1][i] - 9.0l * Phi4[0][i] );
					tmpPhi5( i ) = N5( i ) + dx / 24.0l
									* ( 55.0l * Phi5[3][i] - 59.0l * Phi5[2][i] + 37.0l * Phi5[1][i] - 9.0l * Phi5[0][i] );
				}
				//corrector:
				tmpPhi1 = matrA1 * tmpPhi1;
				tmpPhi2 = matrA1 * tmpPhi2;
				tmpPhi3 = matrA1 * tmpPhi3;
				tmpPhi4 = matrA1 * tmpPhi4;
				tmpPhi5 = matrA1 * tmpPhi5 + vectF1;
				for( int i = 0; i < EQ_NUM; ++i )
				{
					N1( i ) = N1( i ) + dx / 24.0l
									* ( 9.0l * tmpPhi1( i ) + 19.0l * Phi1[3][i] - 5.0l * Phi1[2][i] + Phi1[1][i] );
					N2( i ) = N2( i ) + dx / 24.0l
									* ( 9.0l * tmpPhi2( i ) + 19.0l * Phi2[3][i] - 5.0l * Phi2[2][i] + Phi2[1][i] );
					N3( i ) = N3( i ) + dx / 24.0l
									* ( 9.0l * tmpPhi3( i ) + 19.0l * Phi3[3][i] - 5.0l * Phi3[2][i] + Phi3[1][i] );
					N4( i ) = N4( i ) + dx / 24.0l
									* ( 9.0l * tmpPhi4( i ) + 19.0l * Phi4[3][i] - 5.0l * Phi4[2][i] + Phi4[1][i] );
					N5( i ) = N5( i ) + dx / 24.0l
									* ( 9.0l * tmpPhi5( i ) + 19.0l * Phi5[3][i] - 5.0l * Phi5[2][i] + Phi5[1][i] );
				}
			}
			else
			{
				rungeKutta->calc3( matrA, matrA1, vectF, vectF1, dx, 0, &N1 );
				rungeKutta->calc3( matrA, matrA1, vectF, vectF1, dx, 0, &N2 );
				rungeKutta->calc3( matrA, matrA1, vectF, vectF1, dx, 0, &N3 );
				rungeKutta->calc3( matrA, matrA1, vectF, vectF1, dx, 0, &N4 );
				rungeKutta->calc3( matrA, matrA1, vectF, vectF1, dx, 1, &N5 );
			}

			rgkTime += time( 0 ) - rgkTime1;

			N1orthog = N1;
			N2orthog = N2;
			N3orthog = N3;
			N4orthog = N4;
			N5orthog = N5;

			orthoTime1 = time( 0 );
			orthoBuilder->orthonorm( 1, x, &N1orthog );
			orthoBuilder->orthonorm( 2, x, &N2orthog );
			orthoBuilder->orthonorm( 3, x, &N3orthog );
			orthoBuilder->orthonorm( 4, x, &N4orthog );
			orthoBuilder->orthonorm( 5, x, &N5orthog );
			orthoTime += time( 0 ) - orthoTime1;

			if( orthoBuilder->checkOrtho( x, N2orthog, N3orthog, N4orthog, N5orthog, N2, N3, N4, N5 ) == 1 )	//if the orthonormalization is needed
			{
				active = 1;		//if orthonormalization has been performed, we have to restart the ABM method

				N1 = N1orthog;
				N2 = N2orthog;
				N3 = N3orthog;
				N4 = N4orthog;
				N5 = N5orthog;
				
				orthoBuilder->setOrthoDoneInfo( x );
				cout << " --- at x = " << x << " ortho is needed\n";
			}
			else
			{
				++active;	//if no orthonormalization has been done, we have one more solution that can be used in ABM method
				orthoBuilder->setNextSolVects( x, N1, N2, N3, N4, N5 );
			}
		}

		buildSolnTime1 = time( 0 );
		orthoBuilder->buildSolution( &mesh );
		buildSolnTime += time( 0 ) - buildSolnTime1;

		if( preLin != 0 )
		{
			cont = checkConv();
		}
		++iter;
		//cout << " : " << iter << endl;
	}while( cont == 1 && iter <= MAX_NEWTON_ITER );

	if( iter >= MAX_NEWTON_ITER )
	{
		maxNewtonIterReached = 1;
	}

	for( int x = 0; x < Km; ++x )
	{
		for( int i = 0; i < eq_num; ++i )
		{
			PL_NUM tmpD2N = ( mesh[x].Nk1[i] - mesh[x].Nk0[i] ) / beta / dt / dt - mesh[x].d1N0[i] / beta / dt
						- ( 0.5l - beta ) / beta * mesh[x].d2N0[i];
			PL_NUM tmpD1N = mesh[x].d1N0[i] + 0.5l * dt * ( mesh[x].d2N0[i] + tmpD2N );
			mesh[x].Nk0[i] = mesh[x].Nk1[i];
			mesh[x].d1N0[i] = tmpD1N;
			mesh[x].d2N0[i] = tmpD2N;
		}
	}

	totalTime += time( 0 ) - totalTime1;

	copyToResArr();

	PL_NUM sum = 0.0l;
	for( int y = 0; y < Km; ++y )
	{
		sum += mesh[y].Nk1[1] * mesh[y].Nk1[1];
	}
	//return mesh[ ( Km - 1 ) / 2 ].Nk1[1];
	return sum;
}

template<class PL_NUM>
int Solver<PL_NUM>::checkConv()
{
	//for( int x = 1; x < Km; ++x )
	//{
	//	for( int i = 0; i < eq_num; ++i )
	//	{
	//		if( mesh[x].Nk[i] != 0.0l ) 
	//		{
	//			if( fabs( ( mesh[x].Nk1[i] - mesh[x].Nk[i] ) / mesh[x].Nk[i] ) < ALMOST_ZERO )
	//			{
	//				return 0;
	//			}
	//		}
	//		else
	//		{
	//			if( fabs( mesh[x].Nk1[i] ) < ALMOST_ZERO )
	//			{
	//				return 0;
	//			}
	//		}
	//	}
	//}
	//return 1;

	for( int x = 0; x < Km; ++x )
	{
		for( int i = 0; i < eq_num; ++i )
		{
			if( fabs( mesh[x].Nk[i] ) > ALMOST_ZERO )
			{
				if( fabs( ( mesh[x].Nk1[i] - mesh[x].Nk[i] ) / mesh[x].Nk[i] ) > 0.0000001l )
				{
					return 1;
				}
			}
			else
			{
				if( fabs( mesh[x].Nk1[i] ) > ALMOST_ZERO )
				{
					return 1;
				}
			}
		}
	}
	return 0;
}

template<class PL_NUM>
void Solver<PL_NUM>::dump_sol()
{
	ofstream dumpSol;
	stringstream ss;
	ss << "sol_" << cur_t << ".txt";
	
	dumpSol.open ( ss.str() );
	
	for( int x = 0; x < Km; ++x )
	{
		for( int i = 0; i < eq_num; ++i )
		{
			dumpSol << mesh[x].Nk1[i] << " ";
		}
		dumpSol << endl;
	}


	dumpSol.close();
	return;
}

template<class PL_NUM>
void Solver<PL_NUM>::dump_check_sol( int fNum )
{
	N_PRES t = cur_t;

	int minusOne = -1;
	PL_NUM sum = 0.0;
	for( int i = 0; i <= 100000; ++i )
	{
		PL_NUM omg = (long double)( (long double)M_PI * (long double)M_PI * ( 2.0 * i + 1.0 ) * ( 2.0 * i + 1.0 ) ) * h / 2.0l / a / a * sqrt( B22 / 3.0l / rho );

		minusOne = -minusOne;

		sum = sum + (long double)minusOne / ( 2.0 * i + 1.0 ) / ( 2.0 * i + 1.0 ) / ( 2.0 * i + 1.0 ) / ( 2.0 * i + 1.0 ) / ( 2.0 * i + 1.0 ) * cos( omg * t );
	}
	PL_NUM wTheor;
	wTheor = - p0 * a * a * a * a / h / h / h / B22 * ( 5.0l / 32.0l - 48.0l / M_PI / M_PI / M_PI / M_PI / M_PI * sum );

	stringstream ss;
	if( fNum >= 0 )
	{
		ss << "test_sol_" << fNum << ".txt";
	}
	else
	{
		ss << "test_sol.txt";
	}
	ofstream of1( ss.str(), ofstream::app );
	of1 << cur_t << " ; " << mesh[ ( Km - 1 ) / 2 ].Nk1[1] << " ; " << wTheor /*<< " ; " << resArr[ curTimeStep * Km * eq_num + ( Km - 1 ) / 2 * eq_num + 1]*/
		<< " ; " << fabs( ( wTheor - mesh[ ( Km - 1 ) / 2 ].Nk1[1] ) / wTheor ) << endl;
	of1.close();
}

template<class PL_NUM>
void Solver<PL_NUM>::dump_whole_sol( int var )
{
	stringstream ss;
	ss << "sol_whole_" << curTimeStep << ".bin";
	ofstream of1( ss.str(), ofstream::out | ofstream::binary );
	for( int y = 0; y < Km; ++y )
	{
		of1.write( reinterpret_cast<char*>( &( mesh[y].Nk1[var] ) ), sizeof(PL_NUM) );
	}
	of1.close();
}

template<class PL_NUM>
void Solver<PL_NUM>::dump_left_border_vals()
{
	ofstream of1( "sol_left_border.txt", ofstream::app );
	of1 << "v : " << mesh[0].Nk1[0] << "\tw : " << mesh[0].Nk1[1] << "\tW : " << mesh[0].Nk1[2] << "\tMyy : " << mesh[0].Nk1[5] << "\tNyy : " << mesh[0].Nk1[3] << "\tNyz : " << mesh[0].Nk1[4] << endl;
	of1.close();
}

template<class PL_NUM>
void Solver<PL_NUM>::dumpSolAll( int fNum )
{
	stringstream ss;
	if( fNum >= 0 )
	{
		ss << "all_test_sol_" << fNum << ".txt";
	}
	else
	{
		ss << "all_test_sol.txt";
	}
	ofstream of1( ss.str(), ofstream::app );
	of1 << cur_t;
	for( int i = 0; i < eq_num; ++i )
	{
		of1 << " ; " << mesh[ ( Km - 1 ) / 2 ].Nk1[i];
	}
	of1 << endl;
	of1.close();
}

//template<class PL_NUM>
//void Solver<PL_NUM>::dumpMatrA()
//{
//	ofstream of( "MatrA.txt", ofstream::app );
//
//	for( int i = 0; i < EQ_NUM; ++i )
//	{
//		for( int j = 0; j < EQ_NUM; ++j )
//		{
//			of << nonlin_matr_A( i, j ) << " ";
//		}
//		of << endl;
//	}
//	of << "\n============================================\n";
//}

#endif