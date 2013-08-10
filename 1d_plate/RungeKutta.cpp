#include "RungeKutta.h"

RungeKutta::RungeKutta( int _eq_num )
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

void RungeKutta::calc( const Matrix<PL_NUM, EQ_NUM, EQ_NUM, RowMajor> &A, const Matrix<PL_NUM, EQ_NUM, 1> &f, PL_NUM dx, int hom, Matrix<PL_NUM, EQ_NUM, 1>* x )
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

//void RungeKutta::calc_test( long double *Ni, long double **Ai, long double *fi, long Hom, PL_NUM dm )
//{
//	int i, j;
//	int K = eq_num;
//	int r;
//	long double ym;
//
//	/*
//	long double C1, C2, C3, C4, d21, d32, d31, d43, d42, d41; // parameters in the Runge-Kutta scheme
//	C2 = ( 2.0 * rgk_v - 1.0 ) / 12.0 / rgk_u / ( rgk_v - rgk_u ) / ( 1.0 - rgk_u );
//	C3 = ( 1.0 - 2.0 * rgk_u ) / 12.0 / rgk_v / ( rgk_v - rgk_u ) / ( 1.0 - rgk_v );
//	C4 = ( 6.0 * rgk_u * rgk_v - 4.0 * rgk_u - 4.0 * rgk_v + 3.0 ) / 12.0 / ( 1.0 - rgk_u ) / ( 1.0 - rgk_v );
//	C1 = 1.0 - C2 - C3 - C4;
//	d21 = rgk_u;
//	d32 = 1.0 / 24.0 / C3 / rgk_u / ( 1.0 - rgk_v );
//	d31 = rgk_v - d32;
//	d43 = ( 1.0 - 2.0 * rgk_u ) / 12.0 / C4 / rgk_v / ( rgk_v - rgk_u );
//	d42 = -( rgk_v * ( 4.0 * rgk_v - 5.0 ) - rgk_u + 2.0 ) / 24.0 / C4 / rgk_u / ( rgk_v - rgk_u ) / ( 1.0 - rgk_v );
//	d41 = 1.0 - d42 - d43;*/
//
//	long double *ph1, *ph2, *ph3, *ph4;
//	ph1=new long double[K+2];
//	ph2=new long double[K+2];
//	ph3=new long double[K+2];
//	ph4=new long double[K+2];
//
//	long double *Fi;
//	Fi=new long double[K+2];
//
//	for( r = 1; r <= K; r++ )
//    {
//		ph1[r] = 0.0;
//		ph2[r] = 0.0;
//		ph3[r] = 0.0;
//		ph4[r] = 0.0;
//    }
//	if( Hom == 0 ) // check if system is homogeneous
//	{
//		for( i = 1; i <= K; i++ )
//			for( j = 1; j <= K; j++ )
//				ph1[i] += Ai[i][j] * dm * Ni[j];
//		for( i = 1; i <= K; i++ )
//			for( j = 1; j <= K; j++ )
//				ph2[i] += Ai[i][j] * dm * ( Ni[j] + rgk_d21 * ph1[j] );
//		for( i = 1; i <= K; i++ )
//			for( j = 1; j <= K; j++ )
//				ph3[i] += Ai[i][j] * dm * ( Ni[j] + rgk_d31 * ph1[j] + rgk_d32 * ph2[j] );
//		for( i = 1; i <= K; i++ )
//			for( j = 1; j <= K; j++ )
//				ph4[i] += Ai[i][j] * dm * ( Ni[j] + rgk_d41 * ph1[j] + rgk_d42 * ph2[j] + rgk_d43 * ph3[j] );
//		for( i = 1; i <= K; i++ )
//			Fi[i] = ph1[i];
//
////      for (i=1;i<=K;i++)
////         Ni[i]=Ni[i]+C1*ph1[i]+C2*ph2[i]+C3*ph3[i]+C4*ph4[i];
//	}
//	else
//	{
//
//		for( i = 1; i <= K; i++ )
//			for( j = 1; j <= K; j++ )
//				ph1[i] += Ai[i][j] * dm * Ni[j];
//		for( i = 1; i <= K; i++ )
//			for( j = 1; j <= K; j++ )
//				ph2[i] += Ai[i][j] * dm * ( Ni[j] + rgk_d21 * ph1[j] );
//		for( i = 1; i <= K; i++ )
//			for( j = 1; j <= K; j++ )
//				ph3[i] += Ai[i][j] * dm * ( Ni[j] + rgk_d31 * ph1[j] + rgk_d32 * ph2[j] );
//		for( i = 1; i <= K; i++ )
//			for( j = 1; j <= K; j++ )
//				ph4[i] += Ai[i][j] * dm * ( Ni[j] + rgk_d41 * ph1[j] + rgk_d42 * ph2[j] + rgk_d43 * ph3[j] );
//		for( i = 1; i <= K; i++ )
//		{
//			ph1[i] += fi[i] * dm;
//			ph2[i] += fi[i] * dm;
//			ph3[i] += fi[i] * dm;
//			ph4[i] += fi[i] * dm;
//			Fi[i] = ph1[i];
//		}
///*      for(i=1;i<=K;i++)
//       {
//         ph1[i]=dm*ph1[i];
//         ph2[i]=dm*ph2[i];
//         ph3[i]=dm*ph3[i];
//         ph4[i]=dm*ph4[i];
//         }
//*/
//	}
//
//	for ( i = 1; i <= K; i++ )
//		Ni[i] = Ni[i] + rgk_C1 * ph1[i] + rgk_C2 * ph2[i] + rgk_C3 * ph3[i] + rgk_C4 * ph4[i];
//}
//
//void RungeKutta::do_test()
//{
//	//srand( time( 0 ) );
//
//	PL_NUM dx = 0.000001;
//	int Hom = 1;
//
//	vector<PL_NUM>* A = new vector<PL_NUM>( eq_num * eq_num , 0.0);
//	vector<PL_NUM>* f = new vector<PL_NUM>( eq_num, 0.0);
//	vector<PL_NUM>* x = new vector<PL_NUM>( eq_num, 0.0);
//	vector<PL_NUM>* x1 = new vector<PL_NUM>( eq_num, 0.0);
//
//	long double *Ni;
//	long double **Ai;
//	long double *fi;
//	Ni = new long double [eq_num + 2];
//	fi = new long double [eq_num + 2];
//	Ai = new long double* [eq_num + 2];
//	for( int i = 0; i <= eq_num + 1; i++ )
//	{
//		Ai[i] = new long double [eq_num + 2];
//	}
//
//	for( int i = 0; i < eq_num; ++i ){
//		(*f)[i] = ( PL_NUM )( rand() % 1000 );
//		(*x)[i] = ( PL_NUM )( rand() % 1000 );
//		fi[i + 1] = (*f)[i];
//		Ni[i + 1] = (*x)[i];
//		(*x1)[i] = (*x)[i];
//		for( int j = 0; j < eq_num; ++j ){
//			(*A)[eq_num * i + j] = ( PL_NUM )( rand() % 1000 );
//			Ai[i + 1][j + 1] = (*A)[ eq_num * i + j];
//		}
//	}
//
//	cout << "ddd\n";
//
//	calc( *A, *f, dx, 1, x, 0 );
//	calc_test( Ni, Ai, fi, 1, dx );
//	calc( *A, *f, dx, 1, x1, 0 );
//
//	for( int i = 0; i < eq_num; ++i ){
//		cout << "my\t" << (*f)[i] << "\t\told\t" << fi[i + 1] << endl;
//	}
//	cout << endl;
//	for( int i = 0; i < eq_num; ++i ){
//		cout << "my\t" << (*x)[i] << "\t\told\t" << (*x1)[i] << "\t\told\t" << Ni[i + 1] << endl;
//	}
//}