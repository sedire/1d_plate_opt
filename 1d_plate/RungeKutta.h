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
	~RungeKutta() {};
	inline void calc( const Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor> &A, const Matrix<PL_NUM, EQ_NUM, 1> &f, PL_NUM dx, int hom, Matrix<PL_NUM, EQ_NUM, 1>* x );			//method for solving a system of ODE like dy/dx = Ax + f
	int eq_num;
private:
	//int eq_num;
	PL_NUM rgk_u;
	PL_NUM rgk_v;
	PL_NUM rgk_C1, rgk_C2, rgk_C3, rgk_C4;
	PL_NUM rgk_d21, rgk_d32, rgk_d31, rgk_d43, rgk_d42, rgk_d41;

	Matrix<PL_NUM, EQ_NUM, 1> f1;
	Matrix<PL_NUM, EQ_NUM, 1> f2;
	Matrix<PL_NUM, EQ_NUM, 1> f3;
	Matrix<PL_NUM, EQ_NUM, 1> f4;

	Matrix<PL_NUM, EQ_NUM, 1> Af1;
	Matrix<PL_NUM, EQ_NUM, 1> Af2;
	Matrix<PL_NUM, EQ_NUM, 1> Af3;
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

//template<class PL_NUM>
//void RungeKutta<PL_NUM>::calc( const Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor> &A, const Matrix<PL_NUM, EQ_NUM, 1> &f, PL_NUM dx, int hom, Matrix<PL_NUM, EQ_NUM, 1>* x )
//{
//	//for( int i = 0; i < eq_num; ++i ) {
//	//	f1( i ) = 0.0;
//	//	f2( i ) = 0.0;
//	//	f3( i ) = 0.0;
//	//	f4( i ) = 0.0;
//	//}
//
//	f1 = dx * ( A * (*x) );
//	//for( int i = 0; i < eq_num; ++i ) {				//f1 = dx * Fhi( x )
//	//	for( int j = 0; j < eq_num; ++j ) {
//	//		f1[i] += dx * A[eq_num * i + j] * (x)[j];
//	//	}
//	//}
//
//	Af1 = dx * ( A * f1 );
//	//f2 += dx * ( A * ( (*x) + rgk_d21 * f1 ) );
//	f2 = rgk_d21 * Af1 + f1;
//
//	//for( int i = 0; i < eq_num; ++i ) {				//f2 = dx * Fhi( x + d21 * f1 )
//	//	for( int j = 0; j < eq_num; ++j ) {
//	//		f2[i] += dx * A[eq_num * i + j] * ( (x)[j] + rgk_d21 * f1[j] );
//	//	}
//	//}
//
//	Af2 = dx * ( A * f2 );
//	//f3 += dx * ( A * ( (*x) + rgk_d31 * f1 + rgk_d32 * f2 ) );
//	f3 = rgk_d31 * Af1 + rgk_d32 * Af2 + f1;
//
//	//for( int i = 0; i < eq_num; ++i ) {				//f3 = dx * Fhi( x + d31 * f1 + d32 * f2 )
//	//	for( int j = 0; j < eq_num; ++j ) {
//	//		f3[i] += dx * A[eq_num * i + j] * ( (x)[j] + rgk_d31 * f1[j] + rgk_d32 * f2[j] );
//	//	}
//	//}
//
//	Af3 = dx * ( A * f3 );
//	//f4 += dx * ( A * ( (*x) + rgk_d41 * f1 + rgk_d42 * f2 + rgk_d43 * f3 ) );
//	f4 = rgk_d41 * Af1 + rgk_d42 * Af2 + rgk_d43 * Af3 + f1;
//
//	//for( int i = 0; i < eq_num; ++i ) {				//f2 = dx * Fhi( x + d41 * f1 + d42 * f2 + d43 * f3 )
//	//	for( int j = 0; j < eq_num; ++j ) {
//	//		f4[i] += dx * A[eq_num * i + j] * ( (x)[j] + rgk_d41 * f1[j] + rgk_d42 * f2[j] + rgk_d43 * f3[j] );
//	//	}
//	//}
//	if( hom != 0 ) {
//		f1 += dx * f;
//		f2 += dx * f;
//		f3 += dx * f;
//		f4 += dx * f;
//		//for( int i = 0; i < eq_num; ++i ) {
//		//	f1[i] += dx * f[i];
//		//	f2[i] += dx * f[i];
//		//	f3[i] += dx * f[i];
//		//	f4[i] += dx * f[i];
//		//}
//	}
//
//	//for( int i = 0; i < eq_num; ++i ) {
//	//		(*Phi)[i] = f1[i];
//	//}
//
//	(*x) += rgk_C1 * f1 + rgk_C2 * f2 + rgk_C3 * f3 + rgk_C4 * f4;
//	//for( int i = 0; i < eq_num; ++i ) {
//	//	(x)[i] = (x)[i] + rgk_C1 * f1[i] + rgk_C2 * f2[i] + rgk_C3 * f3[i] + rgk_C4 * f4[i];
//	//}
//}

template<class PL_NUM>
void RungeKutta<PL_NUM>::calc( const Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor> &A, const Matrix<PL_NUM, EQ_NUM, 1> &f, PL_NUM dx, int hom, Matrix<PL_NUM, EQ_NUM, 1>* x )
{
	f1( 0 ) = dx * A( 0, 3 ) * (*x)( 3 );
	f1( 1 ) = dx * A( 1, 2 ) * (*x)( 2 );
	f1( 2 ) = dx * A( 2, 5 ) * (*x)( 5 );
	f1( 3 ) = dx * ( A( 3, 0 ) * (*x)( 0 ) + A( 3, 1 ) * (*x)( 1 ) + A( 3, 2 ) * (*x)( 2 ) + A( 3, 3 ) * (*x)( 3 )
			 + A( 3, 6 ) * (*x)( 6 ) + A( 3, 7 ) * (*x)( 7 ) );
	f1( 4 ) = dx * ( A( 4, 0 ) * (*x)( 0 ) + A( 4, 1 ) * (*x)( 1 ) + A( 4, 2 ) * (*x)( 2 ) + A( 4, 6 ) * (*x)( 6 )
			 + A( 4, 7 ) * (*x)( 7 ) );
	f1( 5 ) = dx * ( A( 5, 1 ) * (*x)( 1 ) + A( 5, 2 ) * (*x)( 2 ) + A( 5, 4 ) * (*x)( 4 ) + A( 5, 5 ) * (*x)( 5 )
			 + A( 5, 6 ) * (*x)( 6 ) + A( 5, 7 ) * (*x)( 7 ) );
	f1( 6 ) = dx * A( 6, 7 ) * (*x)( 7 );
	f1( 7 ) = dx * ( A( 7, 0 ) * (*x)( 0 ) + A( 7, 1 ) * (*x)( 1 ) + A( 7, 6 ) * (*x)( 6 ) + A( 7, 7 ) * (*x)( 7 ) );

	Af1( 0 ) = dx * A( 0, 3 ) * f1( 3 );
	Af1( 1 ) = dx * A( 1, 2 ) * f1( 2 );
	Af1( 2 ) = dx * A( 2, 5 ) * f1( 5 );
	Af1( 3 ) = dx * ( A( 3, 0 ) * f1( 0 ) + A( 3, 1 ) * f1( 1 ) + A( 3, 2 ) * f1( 2 ) + A( 3, 3 ) * f1( 3 )
			 + A( 3, 6 ) * f1( 6 ) + A( 3, 7 ) * f1( 7 ) );
	Af1( 4 ) = dx * ( A( 4, 0 ) * f1( 0 ) + A( 4, 1 ) * f1( 1 ) + A( 4, 2 ) * f1( 2 ) + A( 4, 6 ) * f1( 6 )
			 + A( 4, 7 ) * f1( 7 ) );
	Af1( 5 ) = dx * ( A( 5, 1 ) * f1( 1 ) + A( 5, 2 ) * f1( 2 ) + A( 5, 4 ) * f1( 4 ) + A( 5, 5 ) * f1( 5 )
			 + A( 5, 6 ) * f1( 6 ) + A( 5, 7 ) * f1( 7 ) );
	Af1( 6 ) = dx * A( 6, 7 ) * f1( 7 );
	Af1( 7 ) = dx * ( A( 7, 0 ) * f1( 0 ) + A( 7, 1 ) * f1( 1 ) + A( 7, 6 ) * f1( 6 ) + A( 7, 7 ) * f1( 7 ) );
	f2 = rgk_d21 * Af1 + f1;

	Af2( 0 ) = dx * A( 0, 3 ) * f2( 3 );
	Af2( 1 ) = dx * A( 1, 2 ) * f2( 2 );
	Af2( 2 ) = dx * A( 2, 5 ) * f2( 5 );
	Af2( 3 ) = dx * ( A( 3, 0 ) * f2( 0 ) + A( 3, 1 ) * f2( 1 ) + A( 3, 2 ) * f2( 2 ) + A( 3, 3 ) * f2( 3 )
			 + A( 3, 6 ) * f2( 6 ) + A( 3, 7 ) * f2( 7 ) );
	Af2( 4 ) = dx * ( A( 4, 0 ) * f2( 0 ) + A( 4, 1 ) * f2( 1 ) + A( 4, 2 ) * f2( 2 ) + A( 4, 6 ) * f2( 6 )
			 + A( 4, 7 ) * f2( 7 ) );
	Af2( 5 ) = dx * ( A( 5, 1 ) * f2( 1 ) + A( 5, 2 ) * f2( 2 ) + A( 5, 4 ) * f2( 4 ) + A( 5, 5 ) * f2( 5 )
			 + A( 5, 6 ) * f2( 6 ) + A( 5, 7 ) * f2( 7 ) );
	Af2( 6 ) = dx * A( 6, 7 ) * f2( 7 );
	Af2( 7 ) = dx * ( A( 7, 0 ) * f2( 0 ) + A( 7, 1 ) * f2( 1 ) + A( 7, 6 ) * f2( 6 ) + A( 7, 7 ) * f2( 7 ) );
	f3 = rgk_d31 * Af1 + rgk_d32 * Af2 + f1;


	Af3( 0 ) = dx * A( 0, 3 ) * f3( 3 );
	Af3( 1 ) = dx * A( 1, 2 ) * f3( 2 );
	Af3( 2 ) = dx * A( 2, 5 ) * f3( 5 );
	Af3( 3 ) = dx * ( A( 3, 0 ) * f3( 0 ) + A( 3, 1 ) * f3( 1 ) + A( 3, 2 ) * f3( 2 ) + A( 3, 3 ) * f3( 3 )
			 + A( 3, 6 ) * f3( 6 ) + A( 3, 7 ) * f3( 7 ) );
	Af3( 4 ) = dx * ( A( 4, 0 ) * f3( 0 ) + A( 4, 1 ) * f3( 1 ) + A( 4, 2 ) * f3( 2 ) + A( 4, 6 ) * f3( 6 )
			 + A( 4, 7 ) * f3( 7 ) );
	Af3( 5 ) = dx * ( A( 5, 1 ) * f3( 1 ) + A( 5, 2 ) * f3( 2 ) + A( 5, 4 ) * f3( 4 ) + A( 5, 5 ) * f3( 5 )
			 + A( 5, 6 ) * f3( 6 ) + A( 5, 7 ) * f3( 7 ) );
	Af3( 6 ) = dx * A( 6, 7 ) * f3( 7 );
	Af3( 7 ) = dx * ( A( 7, 0 ) * f3( 0 ) + A( 7, 1 ) * f3( 1 ) + A( 7, 6 ) * f3( 6 ) + A( 7, 7 ) * f3( 7 ) );
	f4 = rgk_d41 * Af1 + rgk_d42 * Af2 + rgk_d43 * Af3 + f1;

	if( hom != 0 ) {
		f1 += dx * f;
		f2 += dx * f;
		f3 += dx * f;
		f4 += dx * f;
	}

	(*x) += rgk_C1 * f1 + rgk_C2 * f2 + rgk_C3 * f3 + rgk_C4 * f4;
}

#endif