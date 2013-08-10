#ifndef _PLATE_1D_RUNGEKUTTA_
#define _PLATE_1D_RUNGEKUTTA_ 1

#include <vector>
#include <iostream>
#include <time.h>
#include "plate_var_types.h"
#include <Eigen/Eigen>
#include <complex>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using namespace Eigen;

class RungeKutta
{
public:
	RungeKutta( int _eq_num );
	~RungeKutta();
	void calc( const Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor> &A, const Matrix<PL_NUM, EQ_NUM, 1> &f, PL_NUM dx, int hom, Matrix<PL_NUM, EQ_NUM, 1>* x );			//method for solving a system of ODE like dy/dx = Ax + f
	//void calc_test( long double *Ni, long double **Ai, long double *fi, long Hom, PL_NUM dm );						//old working one, for testing
	//void do_test();

private:
	int eq_num;
	PL_NUM rgk_u;
	PL_NUM rgk_v;
	PL_NUM rgk_C1, rgk_C2, rgk_C3, rgk_C4;
	PL_NUM rgk_d21, rgk_d32, rgk_d31, rgk_d43, rgk_d42, rgk_d41;

	//vector<PL_NUM> f1;
	//vector<PL_NUM> f2;
	//vector<PL_NUM> f3;
	//vector<PL_NUM> f4;
	//PL_NUM f1[EQ_NUM];
	//PL_NUM f2[EQ_NUM];
	//PL_NUM f3[EQ_NUM];
	//PL_NUM f4[EQ_NUM];
	Matrix<PL_NUM, EQ_NUM, 1> f1;
	Matrix<PL_NUM, EQ_NUM, 1> f2;
	Matrix<PL_NUM, EQ_NUM, 1> f3;
	Matrix<PL_NUM, EQ_NUM, 1> f4;
};

#endif