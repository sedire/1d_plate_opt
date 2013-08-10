#ifndef _PLATE_1D_SOLVER_
#define _PLATE_1D_SOLVER_ 1

#include <time.h>
#include "plate_var_types.h"
#include "Plate.h"
#include "RungeKutta.h"
#include "AdamsBashforth.h"
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

//class SolInfo
//{
//public:
//	SolInfo();
//	~SolInfo();
//
//	vector<PL_NUM> o;		//omega matrix to restore the solution
//	vector<PL_NUM> z1;		//basis vectors of the solution, orthonormalized
//	vector<PL_NUM> z2;
//	vector<PL_NUM> z3;
//	vector<PL_NUM> z4;
//	vector<PL_NUM> z5;
//
//	vector<PL_NUM> C;
//
//	void flushO();
//};

enum {stress_whole, stress_centered};
enum {current_const, current_sin, current_exp_sin};

class Solver
{
public:
	Solver();
	~Solver();

	void loadPlate( Plate* _plate );
	void setTask( PL_NUM _J0, PL_NUM _tauJ );
	void endTask();
	void calcConsts();

	void calc_nonlin_system_run_test( long  _x, long _t );

	PL_NUM do_step();
	void dump_sol();
	void dump_check_sol( int fNum );
	void dump_left_border_vals();
	void dumpMatrA();

	PL_NUM cur_t;
	PL_NUM dt;			//time step
	int curTimeStep;
private:
	int eq_num;				//number of equations in main system

	PL_NUM J0;
	PL_NUM omega;
	PL_NUM  p0;				//constant mechanical load
	int stress_type;
	int current_type;

	PL_NUM B11;
	PL_NUM B22;
	PL_NUM B12;
	PL_NUM By0;
	PL_NUM By1;                                      // in considered boundary-value problem
	PL_NUM By2;   

	PL_NUM eps_0;
	PL_NUM eps_x;
	PL_NUM eps_x_0;

	PL_NUM tauP;
	PL_NUM tauJ;
	PL_NUM rad;

	int Km;				//number of steps by space
	int Kt;				//number of steps by time

	PL_NUM dx;

	PL_NUM al;			//some weird var for normalization. It is said that this will improve the numerical scheme. must be equal to density
	PL_NUM betta;		//parameter at Newmark's time integration scheme

	vector<VarVect> mesh;		//2d mesh for solution.
	//vector<PL_NUM> nonlin_matr_A;		//matrix A for the nonlinear system at certain t and x
	//vector<PL_NUM> nonlin_vect_f;		//vector f on right part of nonlinear system at certain t and x
	//PL_NUM nonlin_matr_A[EQ_NUM * EQ_NUM];		//matrix A for the nonlinear system at certain t and x
	//PL_NUM nonlin_vect_f[EQ_NUM];		//vector f on right part of nonlinear system at certain t and x

	Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor> nonlin_matr_A;
	Matrix<PL_NUM, EQ_NUM, 1> nonlin_vect_f;

	vector<PL_NUM> newmark_A;
	vector<PL_NUM> newmark_B;
	//vector<PL_NUM> lin_matr_A;		//matrix A for the linear system at certain t and x
	//vector<PL_NUM> lin_vect_f;		//vector f on right part of linear system at certain t and x
	PL_NUM lin_matr_A[EQ_NUM * EQ_NUM];		//matrix A for the linear system at certain t and x
	PL_NUM lin_vect_f[EQ_NUM];		//vector f on right part of linear system at certain t and x

	Matrix<PL_NUM, EQ_NUM, 1> N1;
	Matrix<PL_NUM, EQ_NUM, 1> N2;
	Matrix<PL_NUM, EQ_NUM, 1> N3;
	Matrix<PL_NUM, EQ_NUM, 1> N4;
	Matrix<PL_NUM, EQ_NUM, 1> N5;

	//PL_NUM N1[EQ_NUM];//basis vectors of the solution
	//PL_NUM N2[EQ_NUM];
	//PL_NUM N3[EQ_NUM];
	//PL_NUM N4[EQ_NUM];
	//PL_NUM N5[EQ_NUM];


	Plate* plate;
	RungeKutta* rungeKutta;
	AdamsBashforth* adamsBashforth;
	OrthoBuilder* orthoBuilder;

	void calc_Newmark_AB( int _x, int mode );		//don't really know why we need this stuff with mode need to FIX
	void calc_nonlin_system( int _x );
	//void calc_lin_system( int _x );

	//vector<SolInfo> solInfoMap;
	//void orthonorm( int baseV, int n );		//baseV - number of the basis vector to orthonormalize 
											//n - spatial node (i.e, x coordinate). this orthonormalizes N1, ... , N5 and builds Omega
	//void buildSolution();					//builds solution for the current time step
	int checkConv();
};

#endif