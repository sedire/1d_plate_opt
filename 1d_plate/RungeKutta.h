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

template<class PL_NUM>
class RungeKutta
{
public:
	RungeKutta( int _eq_num );
	~RungeKutta();
	void calc( const Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor> &A, const Matrix<PL_NUM, EQ_NUM, 1> &f, PL_NUM dx, int hom, Matrix<PL_NUM, EQ_NUM, 1>* x );			//method for solving a system of ODE like dy/dx = Ax + f

private:
	int eq_num;
	PL_NUM rgk_u;
	PL_NUM rgk_v;
	PL_NUM rgk_C1, rgk_C2, rgk_C3, rgk_C4;
	PL_NUM rgk_d21, rgk_d32, rgk_d31, rgk_d43, rgk_d42, rgk_d41;

	Matrix<PL_NUM, EQ_NUM, 1> f1;
	Matrix<PL_NUM, EQ_NUM, 1> f2;
	Matrix<PL_NUM, EQ_NUM, 1> f3;
	Matrix<PL_NUM, EQ_NUM, 1> f4;
};


template<class PL_NUM>
RungeKutta<PL_NUM>::RungeKutta( int _eq_num )
{
	eq_num = _eq_num;
	rgk_u = 0.3;
	rgk_v = 0.6;

//calculating Runge-Kutta method's parameters
	rgk_C2 = ( 2.0l * rgk_v - 1.0l ) / 12.0l / rgk_u / ( rgk_v - rgk_u ) / ( 1.0l - rgk_u );
	rgk_C3 = ( 1.0l - 2.0l * rgk_u ) / 12.0l / rgk_v / ( rgk_v - rgk_u ) / ( 1.0l - rgk_v );
	rgk_C4 = ( 6.0l * rgk_u * rgk_v - 4.0l * rgk_u - 4.0l * rgk_v + 3.0l ) / 12.0l / ( 1.0l - rgk_u ) / ( 1.0l - rgk_v );
	rgk_C1 = 1.0l - rgk_C2 - rgk_C3 - rgk_C4;
	rgk_d21 = rgk_u;
	rgk_d32 = 1.0l / 24.0l / rgk_C3 / rgk_u / ( 1.0l - rgk_v );
	rgk_d31 = rgk_v - rgk_d32;
	rgk_d43 = ( 1.0l - 2.0l * rgk_u ) / 12.0l / rgk_C4 / rgk_v / ( rgk_v - rgk_u );
	rgk_d42 = -( rgk_v * ( 4.0l * rgk_v - 5.0l ) - rgk_u + 2.0l ) / 24.0l / rgk_C4 / rgk_u / ( rgk_v - rgk_u ) / ( 1.0l - rgk_v );
	rgk_d41 = 1.0l - rgk_d42 - rgk_d43;

	//f1.resize( eq_num );
	//f2.resize( eq_num );
	//f3.resize( eq_num );
	//f4.resize( eq_num );
}

template<class PL_NUM>
void RungeKutta<PL_NUM>::calc( const Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor> &A, const Matrix<PL_NUM, EQ_NUM, 1> &f, PL_NUM dx, int hom, Matrix<PL_NUM, EQ_NUM, 1>* x )
{
	for( int i = 0; i < eq_num; ++i ) {
		f1( i ) = 0.0;
		f2( i ) = 0.0;
		f3( i ) = 0.0;
		f4( i ) = 0.0;
	}

	f1 += dx * ( A * (*x) );
	//for( int i = 0; i < eq_num; ++i ) {				//f1 = dx * Fhi( x )
	//	for( int j = 0; j < eq_num; ++j ) {
	//		f1[i] += dx * A[eq_num * i + j] * (x)[j];
	//	}
	//}
	f2 += dx * ( A * ( (*x) + rgk_d21 * f1 ) );
	//for( int i = 0; i < eq_num; ++i ) {				//f2 = dx * Fhi( x + d21 * f1 )
	//	for( int j = 0; j < eq_num; ++j ) {
	//		f2[i] += dx * A[eq_num * i + j] * ( (x)[j] + rgk_d21 * f1[j] );
	//	}
	//}
	f3 += dx * ( A * ( (*x) + rgk_d31 * f1 + rgk_d32 * f2 ) );
	//for( int i = 0; i < eq_num; ++i ) {				//f3 = dx * Fhi( x + d31 * f1 + d32 * f2 )
	//	for( int j = 0; j < eq_num; ++j ) {
	//		f3[i] += dx * A[eq_num * i + j] * ( (x)[j] + rgk_d31 * f1[j] + rgk_d32 * f2[j] );
	//	}
	//}
	f4 += dx * ( A * ( (*x) + rgk_d41 * f1 + rgk_d42 * f2 + rgk_d43 * f3 ) );
	//for( int i = 0; i < eq_num; ++i ) {				//f2 = dx * Fhi( x + d41 * f1 + d42 * f2 + d43 * f3 )
	//	for( int j = 0; j < eq_num; ++j ) {
	//		f4[i] += dx * A[eq_num * i + j] * ( (x)[j] + rgk_d41 * f1[j] + rgk_d42 * f2[j] + rgk_d43 * f3[j] );
	//	}
	//}
	if( hom != 0 ) {
		f1 += dx * f;
		f2 += dx * f;
		f3 += dx * f;
		f4 += dx * f;
		//for( int i = 0; i < eq_num; ++i ) {
		//	f1[i] += dx * f[i];
		//	f2[i] += dx * f[i];
		//	f3[i] += dx * f[i];
		//	f4[i] += dx * f[i];
		//}
	}

	//for( int i = 0; i < eq_num; ++i ) {
	//		(*Phi)[i] = f1[i];
	//}

	(*x) += rgk_C1 * f1 + rgk_C2 * f2 + rgk_C3 * f3 + rgk_C4 * f4;
	//for( int i = 0; i < eq_num; ++i ) {
	//	(x)[i] = (x)[i] + rgk_C1 * f1[i] + rgk_C2 * f2[i] + rgk_C3 * f3[i] + rgk_C4 * f4[i];
	//}
}

#endif