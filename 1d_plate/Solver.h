#ifndef _PLATE_1D_SOLVER_
#define _PLATE_1D_SOLVER_ 1

#include <time.h>
#include "plate_var_types.h"
#include "Plate.h"
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

using std::cout;
using std::vector;
using std::ofstream;
using std::string;
using std::stringstream;
using std::complex;
using namespace Eigen;

enum {stress_whole, stress_centered};
enum {current_const, current_sin, current_exp_sin};

//template<class PL_NUM>
//struct SolverPar
//{
////parameters of the material
//	PL_NUM E1;				//Young's modulus
//	PL_NUM E2;				//Young's modulus
//	PL_NUM nu21;			//Poisson's ratio	
//	PL_NUM rho;				//composite's density
//
//	PL_NUM sigma_x;			//electric conductivity
//	PL_NUM sigma_x_mu;
//
//	N_PRES h;				//thickness of the plate
//	N_PRES a;				//width of the plate
////other
//	PL_NUM dt;
////stress
//	PL_NUM J0;
//	PL_NUM omega;
//	PL_NUM  p0;				//constant mechanical load
//	int stress_type;
//	int current_type;
//};

template<class PL_NUM>
class Solver
{
public:
	Solver();
	~Solver();

	void setTask( PL_NUM _J0, PL_NUM _tauSin, PL_NUM _tauExp,
				PL_NUM _J0_1, PL_NUM _tauSin_1, PL_NUM _tauExp_1,
				PL_NUM _By0, PL_NUM _p0, PL_NUM _tauP );
	void setMechLoad( PL_NUM _p0, PL_NUM _tauP );

	PL_NUM increaseTime();

	void calc_nonlin_system_run_test( long  _x, long _t );

	PL_NUM do_step();
	void dump_sol();
	void dump_check_sol( int fNum );
	void dump_left_border_vals();
	void dumpMatrA();

	N_PRES cur_t;
	N_PRES dt;			//time step
	int curTimeStep;

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
	int stress_type;
	int current_type;

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
	int Kt;				//number of steps by time

	N_PRES dx;

	N_PRES al;			//some weird var for normalization. It is said that this will improve the numerical scheme. must be equal to density
	PL_NUM betta;		//parameter at Newmark's time integration scheme

	vector<VarVect<PL_NUM> > mesh;		//2d mesh for solution.
	//vector<PL_NUM> nonlin_matr_A;		//matrix A for the nonlinear system at certain t and x
	//vector<PL_NUM> nonlin_vect_f;		//vector f on right part of nonlinear system at certain t and x
	PL_NUM nonlin_matr_A[EQ_NUM][EQ_NUM];		//matrix A for the nonlinear system at certain t and x
	//PL_NUM nonlin_vect_f[EQ_NUM];		//vector f on right part of nonlinear system at certain t and x

	//Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor> nonlin_matr_A;
	Matrix<PL_NUM, EQ_NUM, 1> nonlin_vect_f;

	vector<PL_NUM> newmark_A;
	vector<PL_NUM> newmark_B;

	Matrix<PL_NUM, EQ_NUM, 1> N1;
	Matrix<PL_NUM, EQ_NUM, 1> N2;
	Matrix<PL_NUM, EQ_NUM, 1> N3;
	Matrix<PL_NUM, EQ_NUM, 1> N4;
	Matrix<PL_NUM, EQ_NUM, 1> N5;

	RungeKutta<PL_NUM>* rungeKutta;
	OrthoBuilder<PL_NUM>* orthoBuilder;

	void calc_Newmark_AB( int _x, int mode );		//don't really know why we need this stuff with mode need to FIX
	void calc_nonlin_system( int _x );
	//void calc_lin_system( int _x );

	//vector<SolInfo> solInfoMap;
	//void orthonorm( int baseV, int n );		//baseV - number of the basis vector to orthonormalize 
											//n - spatial node (i.e, x coordinate). this orthonormalizes N1, ... , N5 and builds Omega
	//void buildSolution();					//builds solution for the current time step
	int checkConv();
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
	buildSolnTime( 0 )
{

}

template<class PL_NUM>
Solver<PL_NUM>::~Solver()
{
	delete rungeKutta;
	delete orthoBuilder;
}

template<class PL_NUM>
time_t Solver<PL_NUM>::getOrthoBTime()
{
	return orthoBuilder->orthoTotal;
}

template<class PL_NUM>
void Solver<PL_NUM>::setMechLoad( PL_NUM _p0, PL_NUM _tauP )
{
	tauP = _tauP;//0.01;
	p0 = _p0;//10000000;
	rad = 0.0021 / 100.0;
}

template<class PL_NUM>
void Solver<PL_NUM>::setTask( PL_NUM _J0, PL_NUM _tauSin, PL_NUM _tauExp,
							PL_NUM _J0_1, PL_NUM _tauSin_1, PL_NUM _tauExp_1,
							PL_NUM _By0, PL_NUM _p0, PL_NUM _tauP )
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
	rad = 0.0021 / 100.0;

	_tauSin != 0.0l ? tauSin = _tauSin : tauSin = 1;
	_tauExp != 0.0l ? tauExp = _tauExp : tauExp = 1;
	_tauSin_1 != 0.0l ? tauSin_1 = _tauSin_1 : tauSin_1 = 1;
	_tauExp_1 != 0.0l ? tauExp_1 = _tauExp_1 : tauExp_1 = 1;

	eq_num = 8;
	cur_t = 0.0;
	curTimeStep = 0;

	J0 = _J0;
	J0 *= J0_SCALE;
	J0_1 = _J0_1;
	J0_1 *= J0_SCALE;

	omega = (long double)M_PI / tauSin;
	omega_1 = (long double)M_PI / tauSin_1;
	
	stress_type = stress_centered;
	current_type = current_exp_sin;

	By0 = _By0;
	By0 *= BY0_SCALE;

	eps_0 = 0.000000000008854;
	eps_x = 0.0000000002501502912;

	Km = 10001;
	Kt = 3;

	dt = 0.0001;
	dx = al * a / ( Km - 1 );	//HERE WAS A MISTAKE!!!

	++Km;

	betta = 0.25;

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
			nonlin_matr_A[i][j] = 0.0;
		}
		nonlin_vect_f( i ) = 0.0;
	}

	newmark_A.resize( eq_num, 0.0 );
	newmark_B.resize( eq_num, 0.0 );

	for( int i = 0; i < EQ_NUM; ++i )
	{
		newmark_A[i] = 0.0;
		newmark_B[i] = 0.0;
		N1( i ) = 0.0;
		N2( i ) = 0.0;
		N3( i ) = 0.0;
		N4( i ) = 0.0;
		N5( i ) = 0.0;
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
	By2 = 0.0;
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
void Solver<PL_NUM>::calc_Newmark_AB( int _x, int mode )
{
	if( mode == 0 )
	{
		for( int i = 0; i < eq_num; ++i)
		{
			newmark_A[i] = -mesh[_x].Nk0[i] / betta / dt / dt - mesh[_x].d1N0[i] / betta / dt
							- ( 0.5l - betta ) / betta * mesh[_x].d2N0[i];
			newmark_B[i] = -0.5l * mesh[_x].Nk0[i] / betta / dt + ( 1.0l - 0.5l / betta ) * mesh[_x].d1N0[i]
							- 0.5l * dt * ( 1.0l- 1.0l / betta * ( 0.5l - betta ) ) * mesh[_x].d2N0[i];
		}
		/*newmark_A[0] = -mesh[_x].Nk0[0] / betta / dt / dt - mesh[_x].d1N0[0] / betta / dt
							- ( 0.5 - betta ) / betta * mesh[_x].d2N0[0];
		newmark_A[1] = -mesh[_x].Nk0[1] / betta / dt - mesh[_x].d1N0[1] / betta
							- ( 0.5 - betta ) / betta * dt * mesh[_x].d2N0[1];
		newmark_A[2] = -mesh[_x].Nk0[2] / betta / dt / dt - mesh[_x].d1N0[2] / betta / dt
							- ( 0.5 - betta ) / betta * mesh[_x].d2N0[2];
		newmark_A[7] = -mesh[_x].Nk0[7] / betta / dt / dt - mesh[_x].d1N0[7] / betta / dt
							- ( 0.5 - betta ) / betta * mesh[_x].d2N0[7];

		newmark_B[0] = -0.5 * mesh[_x].Nk0[0] / betta / dt + ( 1.0 - 0.5 / betta ) * mesh[_x].d1N0[0]
							- 0.5 * dt * ( 1.0 - 1.0 / betta * ( 0.5 - betta ) ) * mesh[_x].d2N0[0];
		newmark_B[1] = -0.5 * mesh[_x].Nk0[1] / betta / dt + ( 1.0 - 0.5 / betta ) * mesh[_x].d1N0[1]
							- 0.5 * dt * ( 1.0 - 1.0 / betta * ( 0.5 - betta ) ) * mesh[_x].d2N0[1];
		newmark_B[2] = -0.5 * mesh[_x].Nk0[2] / betta / dt + ( 1.0 - 0.5 / betta ) * mesh[_x].d1N0[2]
							- 0.5 * dt * ( 1.0 - 1.0 / betta * ( 0.5 - betta ) ) * mesh[_x].d2N0[2];
		newmark_B[7] = -0.5 * mesh[_x].Nk0[7] / betta / dt + ( 1.0 - 0.5 / betta ) * mesh[_x].d1N0[7]
							- 0.5 * dt * ( 1.0 - 1.0 / betta * ( 0.5 - betta ) ) * mesh[_x].d2N0[7];*/
	}
	else
	{
		for( int i = 0; i < eq_num; ++i)
		{
			newmark_A[i] = -mesh[_x].Nk1[i] / betta / dt / dt - mesh[_x].d1N[i] / betta / dt
							- ( 0.5l - betta ) / betta * mesh[_x].d2N[i];
			newmark_B[i] = -0.5l * mesh[_x].Nk1[i] / betta / dt + ( 1.0l - 0.5l / betta ) * mesh[_x].d1N[i]
							- 0.5l * dt * ( 1.0l - 1.0l / betta * ( 0.5l - betta ) ) * mesh[_x].d2N[i];
		}
		/*newmark_A[0] = -mesh[_x].Nk1[0] / betta / dt / dt - mesh[_x].d1N[0] / betta / dt
							- ( 0.5 - betta ) / betta * mesh[_x].d2N[0];
		newmark_A[1] = -mesh[_x].Nk1[1] / betta / dt - mesh[_x].d1N[1] / betta
							- ( 0.5 - betta ) / betta * dt * mesh[_x].d2N[1];
		newmark_A[2] = -mesh[_x].Nk1[2] / betta / dt / dt - mesh[_x].d1N[2] / betta / dt
							- ( 0.5 - betta ) / betta * mesh[_x].d2N[2];
		newmark_A[7] = -mesh[_x].Nk1[7] / betta / dt / dt - mesh[_x].d1N[7] / betta / dt
							- ( 0.5 - betta ) / betta * mesh[_x].d2N[7];

		newmark_B[0] = -0.5 * mesh[_x].Nk1[0] / betta / dt + ( 1.0 - 0.5 / betta ) * mesh[_x].d1N[0]
							- 0.5 * dt * ( 1.0 - 1.0 / betta * ( 0.5 - betta ) ) * mesh[_x].d2N[0];
		newmark_B[1] = -0.5 * mesh[_x].Nk1[1] / betta / dt + ( 1.0 - 0.5 / betta ) * mesh[_x].d1N[1]
							- 0.5 * dt * ( 1.0 - 1.0 / betta * ( 0.5 - betta ) ) * mesh[_x].d2N[1];
		newmark_B[2] = -0.5 * mesh[_x].Nk1[2] / betta / dt + ( 1.0 - 0.5 / betta ) * mesh[_x].d1N[2]
							- 0.5 * dt * ( 1.0 - 1.0 / betta * ( 0.5 - betta ) ) * mesh[_x].d2N[2];
		newmark_B[7] = -0.5 * mesh[_x].Nk1[7] / betta / dt + ( 1.0 - 0.5 / betta ) * mesh[_x].d1N[7]
							- 0.5 * dt * ( 1.0 - 1.0 / betta * ( 0.5 - betta ) ) * mesh[_x].d2N[7];*/
	}
}

template<class PL_NUM>
void Solver<PL_NUM>::calc_nonlin_system( int _x )
{
	PL_NUM Jx = 0.0;
	if( current_type == current_const )
	{
		Jx = J0;
	}
	else if( current_type == current_sin )
	{
		Jx = J0 * sin( omega * ( cur_t + dt ) );
	}
	else if( current_type == current_exp_sin )
	{
		if( cur_t + dt < SWITCH_TIME )
		{
			Jx = J0 * exp( -( cur_t + dt ) / tauExp ) * sin( omega * ( cur_t + dt ) );
		}
		else
		{
			Jx = J0_1 * exp( -( cur_t + dt ) / tauExp_1 ) * sin( omega_1 * ( cur_t + dt ) );
		}
	}
	PL_NUM Pimp = 0.0l;
	if( stress_type == stress_centered )
	{
		if( cur_t + dt < tauP && fabs( (long double)_x * dx - a / 2.0 ) < rad )
		{
			Pimp = p0 * sqrt( 1.0l - fabs( (long double)_x * dx - a / 2.0l ) * fabs( (long double)_x * dx - a / 2.0 ) / rad / rad	) 
				* sin( (long double)M_PI * ( cur_t + dt ) / tauP );
		}
	}
	else if( stress_type == stress_whole )
	{
		Pimp = p0;
	}

	//long r;

	nonlin_matr_A[ 0][ 3] = 1.0l / al / h / B22;
	nonlin_matr_A[ 1][ 2] = 1.0l / al;
	nonlin_matr_A[ 2][5 ] = -12.0l /al / h / h / h / B22;

	nonlin_matr_A[ 3][0 ] = rho / al * h / betta / dt / dt + sigma_x * h / 2.0l / betta / dt / al * mesh[_x].Nk[7] * mesh[_x].Nk[7];
	nonlin_matr_A[ 3][1 ] = -sigma_x * h / 4.0l / betta / dt / al * By1 * mesh[_x].Nk[7];
	nonlin_matr_A[ 3][2 ] = -eps_x_0 / 4.0l / betta / dt / al * h * By1 * mesh[_x].Nk[6];
	nonlin_matr_A[ 3][3 ] = eps_x_0 / 2.0l / betta / dt / B22 / al * h * mesh[_x].Nk[6] * mesh[_x].Nk[7];
	nonlin_matr_A[ 3][6 ] = 1.0l / al * ( sigma_x * h * mesh[_x].Nk[7] + eps_x_0 * h / B22 * mesh[_x].Nk[7] * ( 1.0 / 2.0l / betta / dt 
									* mesh[_x].Nk[3]
									+ newmark_B[3] ) - eps_x_0 / 2.0l * h * By1 * ( 1.0l / 2.0l / betta / dt * mesh[_x].Nk[2] + newmark_B[2] ) );
	nonlin_matr_A[ 3][7 ] = 1.0l / al * ( sigma_x * h * mesh[_x].Nk[6] + 2.0l * sigma_x * h * mesh[_x].Nk[7] 
									* ( 1.0 / 2.0l / betta / dt * mesh[_x].Nk[0] + newmark_B[0] )
									- sigma_x * h / 2.0l * By1 * ( 1.0l / 2.0l / betta / dt * mesh[_x].Nk[1] + newmark_B[1] ) 
									+ eps_x_0 * h / B22 * mesh[_x].Nk[6]
									* ( 1.0 / 2.0l / betta / dt * mesh[_x].Nk[3] 
									+ newmark_B[3] ) + h * Jx );

	nonlin_matr_A[ 4][0 ] = -sigma_x * h / 4.0l / betta / dt / al * By1 * mesh[_x].Nk[7];
	nonlin_matr_A[ 4][1 ] = rho / al * h / betta / dt / dt + sigma_x * h / 8.0l / betta / al / dt * ( By1 * By1 + 1.0 / 3.0l * By2 * By2 );
	nonlin_matr_A[ 4][2 ] = 1.0 / 2.0l / betta / al * ( sigma_x * h * h / 12.0l * By2 * mesh[_x].Nk[7] 
									- eps_x_0 * h * mesh[_x].Nk[6] * mesh[_x].Nk[7] ) / dt;
	nonlin_matr_A[ 4][6 ] = -1.0l / al * ( sigma_x * h / 2.0l * By1 + eps_x_0 * h 
									* mesh[_x].Nk[7] * ( 1.0 / 2.0l / betta * mesh[_x].Nk[2] / dt + newmark_B[2] ) );
	nonlin_matr_A[ 4][7 ] = -1.0l / al * ( sigma_x * h / 2.0l * By1 * ( 1.0 / 2.0l / betta * mesh[_x].Nk[0] / dt 
									+ newmark_B[0] ) - ( sigma_x / 12.0l * h * h * By2
									- eps_x_0 * h * mesh[_x].Nk[6] ) * ( 1.0 / 2.0l / betta * mesh[_x].Nk[2] / dt + newmark_B[2] ) );

	nonlin_matr_A[ 5][1 ] = -sigma_x * h * h / 24.0l / betta / dt / al * By2 * mesh[_x].Nk[7];
	nonlin_matr_A[ 5][2 ] = -1.0 / 2.0l / betta / dt / al * ( sigma_x * h * h * h / 12.0l * mesh[_x].Nk[7] 
									* mesh[_x].Nk[7] + eps_x_0 / 12.0l * h * h * By2 * mesh[_x].Nk[6] ) - h * h * h / 12.0l / betta / dt / dt * rho / al;
	nonlin_matr_A[ 5][4 ] = 1.0l / al;
	nonlin_matr_A[ 5][5 ] = eps_x_0 / 2.0l / betta / B22 / dt / al * mesh[_x].Nk[6] * mesh[_x].Nk[7];
	nonlin_matr_A[ 5][6 ] = -eps_x_0 / al * ( h * h / 12.0l * By2 * ( 1.0l / 2.0l / betta / dt * mesh[_x].Nk[2] + newmark_B[2] )
									- mesh[_x].Nk[7] / B22 * ( 1.0l / 2.0l / betta / dt * mesh[_x].Nk[5] + newmark_B[5] ) );
	nonlin_matr_A[ 5][7 ] = -1.0l / al * ( sigma_x * h * h / 12.0l * By2 * ( 1.0 / 2.0l / betta / dt * mesh[_x].Nk[1] + newmark_B[1])
									+ sigma_x * h * h * h / 6.0l * mesh[_x].Nk[7] * ( 1.0l / 2.0l / betta / dt * mesh[_x].Nk[2] 
									+ newmark_B[2] ) - eps_x_0 / B22 * mesh[_x].Nk[6] * ( 1.0l / 2.0l / betta / dt * mesh[_x].Nk[5] + newmark_B[5] ) );

	nonlin_matr_A[ 6][7 ] = 1.0l / 2.0l / ( betta * al * dt );

	nonlin_matr_A[ 7][0 ] = sigma_x_mu / 2.0l / betta / ( dt * al ) * mesh[_x].Nk[7];
	nonlin_matr_A[ 7][1 ] = -sigma_x_mu / 4.0l / betta / ( dt * al ) * By1;
	nonlin_matr_A[ 7][6 ] = sigma_x_mu / al;
	nonlin_matr_A[ 7][7 ] = sigma_x_mu / al * ( 1.0 / 2.0l / betta / dt * mesh[_x].Nk[0] + newmark_B[0]);

	nonlin_vect_f( 3 ) = rho / al * h * newmark_A[0] + 1.0l / al * ( -sigma_x * h * mesh[_x].Nk[6] 
						* mesh[_x].Nk[7] - sigma_x * h / betta / dt * mesh[_x].Nk[7] 
						* mesh[_x].Nk[7] * mesh[_x].Nk[0] - sigma_x * h * mesh[_x].Nk[7] 
						* mesh[_x].Nk[7] * newmark_B[0] 
						+ sigma_x * h / 4.0l / betta / dt * By1 * mesh[_x].Nk[7] * mesh[_x].Nk[1] 
						- eps_x_0 * h / betta / dt / B22 * mesh[_x].Nk[6] * mesh[_x].Nk[7] * mesh[_x].Nk[3] - eps_x_0 * h / B22 
						* mesh[_x].Nk[6] * mesh[_x].Nk[7] * newmark_B[3] 
						+ eps_x_0 * h / 4.0l / betta / dt * By1 * mesh[_x].Nk[6] * mesh[_x].Nk[2] );
	nonlin_vect_f( 4 ) = rho / al * h * newmark_A[1] + Pimp / al + 1.0l / al * ( sigma_x * h / 4.0l / betta * By1 / dt	//55454 h
						* mesh[_x].Nk[7] * mesh[_x].Nk[0] + sigma_x * h / 4.0l * ( By1 * By1 + 1.0l / 3.0l * By2 * By2 ) * newmark_B[1] 
						- sigma_x * h * h / 24.0l / betta / dt * By2 * mesh[_x].Nk[7] * mesh[_x].Nk[2] + eps_x_0 * h / betta / dt 
						* mesh[_x].Nk[6] * mesh[_x].Nk[7] * mesh[_x].Nk[2] + eps_x_0 * h * mesh[_x].Nk[6] 
						* mesh[_x].Nk[7] * newmark_B[2] - h / 2.0l * Jx * By1 );			//By1
	nonlin_vect_f( 5 ) = 1.0l / al * ( sigma_x * h * h / 24.0l / betta / dt * By2 * mesh[_x].Nk[7] * mesh[_x].Nk[1] 
						+ sigma_x * h * h * h / 12.0l / betta / dt * mesh[_x].Nk[7] * mesh[_x].Nk[7] * mesh[_x].Nk[2] 
						+ sigma_x * h * h * h / 12.0l * mesh[_x].Nk[7] * mesh[_x].Nk[7] * newmark_B[2] + eps_x_0 / 24.0l / betta / dt * h * h 
						* By2 * mesh[_x].Nk[6] * mesh[_x].Nk[2] - eps_x_0 / betta / dt / B22 * mesh[_x].Nk[6] 
						* mesh[_x].Nk[7] * mesh[_x].Nk[5] - eps_x_0 / B22 * mesh[_x].Nk[6] * mesh[_x].Nk[7] 
						* newmark_B[5] ) - h * h * h / 12.0l * newmark_A[2] * rho  / al;
	nonlin_vect_f( 6 ) = newmark_B[7] / al;
	nonlin_vect_f( 7 ) = 1.0l / al * ( By2 / h - sigma_x_mu / 2.0l / betta / dt * mesh[_x].Nk[7] * mesh[_x].Nk[0] - 0.5l * sigma_x_mu * By1 * newmark_B[1] );
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
		if( iter == 0 && curTimeStep == 0 )
		{
			preLin = 0;
		}
		else
		{
			preLin = 1;
		}
		calc_Newmark_AB( 0, 1 );

		//solInfoMap[0].flushO();
		orthoBuilder->flushO( 0 );

		N1( 0 ) = 0.0; N1( 1 ) = 0.0; N1( 2 ) = 1.0; N1( 3 ) = 0.0; N1( 4 ) = 0.0; N1( 5 ) = 0.0; N1( 6 ) = 0.0; N1( 7 ) = 0.0;
		N2( 0 ) = 0.0; N2( 1 ) = 0.0; N2( 2 ) = 0.0; N2( 3 ) = 1.0; N2( 4 ) = 0.0; N2( 5 ) = 0.0; N2( 6 ) = 0.0; N2( 7 ) = 0.0;
		N3( 0 ) = 0.0; N3( 1 ) = 0.0; N3( 2 ) = 0.0; N3( 3 ) = 0.0; N3( 4 ) = 1.0; N3( 5 ) = 0.0; N3( 6 ) = 0.0; N3( 7 ) = 0.0;
		N4( 0 ) = 0.0; N4( 1 ) = 0.0; N4( 2 ) = 0.0; N4( 3 ) = 0.0; N4( 4 ) = 0.0; N4( 5 ) = 0.0; 
		if( preLin == 0 )
		{
			N4( 6 ) = 0.0;
		}
		else
		{
			N4( 6 ) = mesh[0].Nk1[0] / 2.0l / betta / dt + newmark_B[0];
		}
		N4( 7 ) = -1.0;
		N5( 0 ) = 0.0; N5( 1 ) = 0.0; N5( 2 ) = 0.0; N5( 3 ) = 0.0; N5( 4 ) = 0.0; N5( 5 ) = 0.0; 
		if( preLin == 0 )
		{
			N5( 6 ) = 0.0;
			N5( 7 ) = 0.0;
		}
		else
		{
			N5( 6 ) = ( newmark_B[0] * mesh[0].Nk1[7] - newmark_B[1] * By0 );
			N5( 7 ) = 0.0l - mesh[0].Nk1[7];
		}

		/*for( int i = 0; i < eq_num; ++i )
		{
			solInfoMap[0].z1[i] = N1[i];
			solInfoMap[0].z2[i] = N2[i];
			solInfoMap[0].z3[i] = N3[i];
			solInfoMap[0].z4[i] = N4[i];
			solInfoMap[0].z5[i] = N5[i];
		}*/
		orthoBuilder->setInitVects( N1, N2, N3, N4, N5 );

		for( int x = 0; x < Km; ++x )
		{
			for( int i = 0; i < eq_num; ++i )
			{
				mesh[x].Nk[i] = mesh[x].Nk1[i];
			}
		}

		for( int x = 0; x < Km - 1; ++x )
		{
			orthoBuilder->flushO( x + 1 );

			calc_Newmark_AB( x, 0 );
			matrTime1 = time( 0 );
			calc_nonlin_system( x );
			matrTime += time( 0 ) - matrTime1;

			rgkTime1 = time( 0 );

			//cout << " x is " << x << endl;

			rungeKutta->calc( nonlin_matr_A, nonlin_vect_f, dx, 0, &N1 );
			rungeKutta->calc( nonlin_matr_A, nonlin_vect_f, dx, 0, &N2 );
			rungeKutta->calc( nonlin_matr_A, nonlin_vect_f, dx, 0, &N3 );
			rungeKutta->calc( nonlin_matr_A, nonlin_vect_f, dx, 0, &N4 );
			rungeKutta->calc( nonlin_matr_A, nonlin_vect_f, dx, 1, &N5 );

			rgkTime += time( 0 ) - rgkTime1;

			orthoTime1 = time( 0 );
			orthoBuilder->orthonorm( 1, x, &N1 );
			orthoBuilder->orthonorm( 2, x, &N2 );
			orthoBuilder->orthonorm( 3, x, &N3 );
			orthoBuilder->orthonorm( 4, x, &N4 );
			orthoBuilder->orthonorm( 5, x, &N5 );
			orthoTime += time( 0 ) - orthoTime1;
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
	}while( cont == 1 );

	//cout << "approximation to the solution on the time step done in " << iter << " steps\n";

	for( int x = 0; x < Km; ++x )
	{
		for( int i = 0; i < eq_num; ++i )
		{
			mesh[x].d2N[i] = ( mesh[x].Nk1[i] - mesh[x].Nk0[i] ) / betta / dt / dt - mesh[x].d1N0[i] / betta / dt
						- ( 0.5l - betta ) / betta * mesh[x].d2N0[i];
			mesh[x].d1N[i] = mesh[x].d1N0[i] + 0.5l * dt * ( mesh[x].d2N0[i] + mesh[x].d2N[i] );
			mesh[x].Nk0[i] = mesh[x].Nk1[i];
			mesh[x].d1N0[i] = mesh[x].d1N[i];
			mesh[x].d2N0[i] = mesh[x].d2N[i];
		}
	}

	totalTime += time( 0 ) - totalTime1;

	return mesh[ ( Km - 1 ) / 2 ].Nk1[1];
}

template<class PL_NUM>
int Solver<PL_NUM>::checkConv()
{
	for( int x = 1; x < Km; ++x )
	{
		for( int i = 0; i < eq_num; ++i )
		{
			if( mesh[x].Nk[i] != 0.0l ) 
			{
				if( fabs( ( mesh[x].Nk1[i] - mesh[x].Nk[i] ) / mesh[x].Nk[i] ) < ALMOST_ZERO )
				{
					return 0;
				}
			}
			else
			{
				if( fabs( mesh[x].Nk1[i] ) < ALMOST_ZERO )
				{
					return 0;
				}
			}
		}
	}
	return 1;
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

	/*int minusOne = -1;
	PL_NUM sum = 0.0;
	for( int i = 0; i <= 1000000; ++i )
	{
		PL_NUM omg = (long double)( (long double)M_PI * (long double)M_PI * ( 2 * i + 1 ) * ( 2 * i + 1 ) ) * h / 2.0l / a / a * sqrt( B22 / 3.0l / rho );

		minusOne = -minusOne;

		sum = sum + (long double)minusOne / ( 2 * i + 1 ) / ( 2 * i + 1 ) / ( 2 * i + 1 ) / ( 2 * i + 1 ) / ( 2 * i + 1 ) * cos( omg * t );
	}
	PL_NUM wTheor;
	wTheor = - p0 * a * a * a * a / h / h / h / B22 * ( 5.0l / 32.0l - 48.0l / M_PI / M_PI / M_PI / M_PI / M_PI * sum );*/

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
	of1 << cur_t << " ; " << mesh[ ( Km - 1 ) / 2 ].Nk1[1] <</* " ; " << wTheor << " ; " << fabs( ( wTheor - mesh[ ( Km - 1 ) / 2 ].Nk1[1] ) / wTheor ) << */endl;
	of1.close();
}

template<class PL_NUM>
void Solver<PL_NUM>::dump_left_border_vals()
{
	ofstream of1( "sol_left_border.txt", ofstream::app );
	of1 << "v : " << mesh[0].Nk1[0] << "\tw : " << mesh[0].Nk1[1] << "\tW : " << mesh[0].Nk1[2] << "\tMyy : " << mesh[0].Nk1[5] << "\tNyy : " << mesh[0].Nk1[3] << "\tNyz : " << mesh[0].Nk1[4] << endl;
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