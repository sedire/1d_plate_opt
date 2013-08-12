#ifndef _PLATE_1D_ORTHOBUILDER_
#define _PLATE_1D_ORTHOBUILDER_ 1

#include "plate_var_types.h"
#include "VarVect.h"
#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <complex>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using namespace Eigen;

template<class PL_NUM>
class SolInfo
{
public:
	SolInfo();
	~SolInfo();

	vector<PL_NUM> o;		//omega matrix to restore the solution
	vector<PL_NUM> z1;		//basis vectors of the solution, orthonormalized
	vector<PL_NUM> z2;
	vector<PL_NUM> z3;
	vector<PL_NUM> z4;
	vector<PL_NUM> z5;

	vector<PL_NUM> C;

	void flushO();
};

template<class PL_NUM>
class OrthoBuilder
{
public:
	OrthoBuilder();
	virtual ~OrthoBuilder();
	virtual void setParams( int _Km );
	//virtual void orthonorm( int baseV, int n, vector<PL_NUM>* NtoOrt ) {};
	virtual void orthonorm( int baseV, int n, Matrix<PL_NUM, EQ_NUM, 1>* NtoOrt ) {};
	virtual void buildSolution( vector<VarVect<PL_NUM> >* _mesh ) {};
	virtual void flushO( int x );
	virtual void setInitVects( const Matrix<PL_NUM, EQ_NUM, 1> &N1, const Matrix<PL_NUM, EQ_NUM, 1> &N2, const Matrix<PL_NUM, EQ_NUM, 1> &N3, const Matrix<PL_NUM, EQ_NUM, 1> &N4, const Matrix<PL_NUM, EQ_NUM, 1> &N5 );
	virtual void orthogTest( int x );
protected:
	int eq_num;
	int Km;
	vector<SolInfo<PL_NUM> > solInfoMap;
};

template<class PL_NUM>
class OrthoBuilderGSh : public OrthoBuilder<PL_NUM>
{
public:
	OrthoBuilderGSh() {};
	~OrthoBuilderGSh() {};
	//void orthonorm( int baseV, int n, vector<PL_NUM>* NtoOrt );
	void orthonorm( int baseV, int n, Matrix<PL_NUM, EQ_NUM, 1>* NtoOrt );
	void buildSolution( vector<VarVect<PL_NUM> >* _mesh );
};


template<class PL_NUM>
SolInfo<PL_NUM>::SolInfo()
{
	o.resize( EQ_NUM * EQ_NUM, 0.0 );
	z1.resize( EQ_NUM );
	z2.resize( EQ_NUM );
	z3.resize( EQ_NUM );
	z4.resize( EQ_NUM );
	z5.resize( EQ_NUM );

	C.resize( EQ_NUM / 2 );
}

template<class PL_NUM>
void SolInfo<PL_NUM>::flushO()
{
	for( int i = 0; i < o.size(); ++i )
	{
		o[i] = 0.0;
	}
}

template<class PL_NUM>
SolInfo<PL_NUM>::~SolInfo()
{

}

template<class PL_NUM>
OrthoBuilder<PL_NUM>::OrthoBuilder()
{
	eq_num = EQ_NUM;
}

template<class PL_NUM>
OrthoBuilder<PL_NUM>::~OrthoBuilder()
{

}

template<class PL_NUM>
void OrthoBuilder<PL_NUM>::flushO( int x )
{
	solInfoMap[x].flushO();
}

template<class PL_NUM>
void OrthoBuilder<PL_NUM>::setParams( int _Km )
{
	Km = _Km;

	solInfoMap.resize( Km );
}

template<class PL_NUM>
void OrthoBuilder<PL_NUM>::setInitVects( const Matrix<PL_NUM, EQ_NUM, 1> &N1, const Matrix<PL_NUM, EQ_NUM, 1> &N2, const Matrix<PL_NUM, EQ_NUM, 1> &N3, const Matrix<PL_NUM, EQ_NUM, 1> &N4, const Matrix<PL_NUM, EQ_NUM, 1> &N5 )
{
	for( int i = 0; i < eq_num; ++i )
	{
		solInfoMap[0].z1[i] = N1( i );
		solInfoMap[0].z2[i] = N2( i );
		solInfoMap[0].z3[i] = N3( i );
		solInfoMap[0].z4[i] = N4( i );
		solInfoMap[0].z5[i] = N5( i );
	}
}

template<class PL_NUM>
void OrthoBuilder<PL_NUM>::orthogTest( int x )
{
	PL_NUM sq1 = 0;
	PL_NUM sq2 = 0;
	PL_NUM sq3 = 0;
	PL_NUM sq4 = 0;
	PL_NUM sq5 = 0;

	PL_NUM sum1 = 0;
	PL_NUM sum2 = 0;
	PL_NUM sum3 = 0;
	PL_NUM sum4 = 0;
	PL_NUM sum5 = 0;

	for( int i = 0; i < eq_num; ++i )
	{
		sq1 += solInfoMap[x].z1[i] * solInfoMap[x].z1[i];
		sq2 += solInfoMap[x].z2[i] * solInfoMap[x].z2[i];
		sq3 += solInfoMap[x].z3[i] * solInfoMap[x].z3[i];
		sq4 += solInfoMap[x].z4[i] * solInfoMap[x].z4[i];
		sq5 += solInfoMap[x].z5[i] * solInfoMap[x].z5[i];
	}

	for( int i = 0; i < eq_num; ++i )
	{
		sum1 += solInfoMap[x].z1[i] * solInfoMap[x].z2[i];
		sum2 += solInfoMap[x].z2[i] * solInfoMap[x].z3[i];
		sum3 += solInfoMap[x].z3[i] * solInfoMap[x].z4[i];
		sum4 += solInfoMap[x].z4[i] * solInfoMap[x].z5[i];
		sum5 += solInfoMap[x].z5[i] * solInfoMap[x].z1[i];
	}

	//cout << "  " << sqrtl( sq1 ) << "  " << sqrtl( sq2 ) << "  " << sqrtl( sq3 ) << "  " << sqrtl( sq4 ) << "  " << sqrtl( sq5 ) << endl;
	cout << "  " << sum1 << "  " << sum2 << "  " << sum3 << "  " << sum4 << "  " << sum5 << endl;

}

template<class PL_NUM>
void OrthoBuilderGSh<PL_NUM>::orthonorm( int baseV, int n, Matrix<PL_NUM, EQ_NUM, 1>* NtoOrt )
{
	long double k11 = 1.414213562373095;
	PL_NUM norm = 0;
	vector<PL_NUM> omega2;
	omega2.resize( eq_num * eq_num, 0.0 );

	if( baseV < 1 || baseV > 5 || n < 0 || n > Km - 1 )
	{
		cout << "Error in orthonorm: bad input\n";
		return;
	}

	if( baseV == 1 )
	{
		for( int i = 0; i < eq_num; ++i )
		{
			solInfoMap[n + 1].o[0 * eq_num + 0] += (*NtoOrt)( i ) * (*NtoOrt)( i );
		}
		solInfoMap[n + 1].o[0 * eq_num + 0] = sqrt( solInfoMap[n + 1].o[0 * eq_num + 0] );
		for( int i = 0; i < eq_num; ++i )
		{
			(*NtoOrt)( i ) = (*NtoOrt)( i ) / solInfoMap[n + 1].o[0 * eq_num + 0];
			solInfoMap[n + 1].z1[i] = (*NtoOrt)( i );
		}
	}
	else if( baseV == 2 )
	{
		for( int i = 0; i < eq_num; ++i )
		{
			norm += (*NtoOrt)( i ) * (*NtoOrt)( i );
			solInfoMap[n + 1].o[1 * eq_num + 0] += (*NtoOrt)( i ) * solInfoMap[n + 1].z1[i];
		}
		norm = sqrt( norm );
		for( int i = 0; i < eq_num; ++i )
		{
			(*NtoOrt)( i ) = (*NtoOrt)( i ) - solInfoMap[n + 1].o[1 * eq_num + 0] * solInfoMap[n + 1].z1[i];
		}
		for( int i = 0; i < eq_num; ++i )
		{
			solInfoMap[n + 1].o[1 * eq_num + 1] += (*NtoOrt)( i ) * (*NtoOrt)( i );
		}
		solInfoMap[n + 1].o[1 * eq_num + 1] = sqrt( solInfoMap[n + 1].o[1 * eq_num + 1] );

		if( ( norm / solInfoMap[n + 1].o[1 * eq_num + 1] ).real() <= k11 )
		{
			for( int i = 0; i < eq_num; ++i )
			{
				(*NtoOrt)( i ) = (*NtoOrt)( i )	/ solInfoMap[n + 1].o[1 * eq_num + 1];
				solInfoMap[n + 1].z2[i] = (*NtoOrt)( i );
			}
		}
		else
		{
			for( int i = 0; i < eq_num; ++i )
			{
				omega2[ 1 * eq_num + 0 ] += (*NtoOrt)( i ) * solInfoMap[n + 1].z1[i];
			}
			for( int i = 0; i < eq_num; ++i )
			{
				(*NtoOrt)( i ) -= omega2[ 1 * eq_num + 0 ] * solInfoMap[n + 1].z1[i];
			}
			for( int i = 0; i < eq_num; ++i )
			{
				omega2[ 1 * eq_num + 1 ] += (*NtoOrt)( i ) * (*NtoOrt)( i );
			}
			omega2[ 1 * eq_num + 1 ] = sqrt( omega2[ 1 * eq_num + 1 ] );

			for( int i = 0; i < eq_num; ++i )
			{
				solInfoMap[n + 1].z2[i] = (*NtoOrt)( i ) / omega2[ 1 * eq_num + 1 ];
				(*NtoOrt)( i ) = solInfoMap[n + 1].z2[i];
			}
			solInfoMap[n + 1].o[1 * eq_num + 0] += omega2[ 1 * eq_num + 0 ];
			solInfoMap[n + 1].o[1 * eq_num + 1] = omega2[ 1 * eq_num + 1 ];
		}
	}
	else if( baseV == 3 )
	{
		for( int i = 0; i < eq_num; ++i )
		{
			norm += (*NtoOrt)( i ) * (*NtoOrt)( i );
			solInfoMap[n + 1].o[2 * eq_num + 0] += (*NtoOrt)( i ) * solInfoMap[n + 1].z1[i];
		}
		norm = sqrt( norm );
		for( int i = 0; i < eq_num; ++i )
		{
			(*NtoOrt)( i ) = (*NtoOrt)( i ) - solInfoMap[n + 1].o[2 * eq_num + 0] * solInfoMap[n + 1].z1[i];
		}
		for( int i = 0; i < eq_num; ++i )
		{
			solInfoMap[n + 1].o[2 * eq_num + 1] += (*NtoOrt)( i ) * solInfoMap[n + 1].z2[i];
		}
		for( int i = 0; i < eq_num; ++i )
		{
			(*NtoOrt)( i ) = (*NtoOrt)( i ) - solInfoMap[n + 1].o[2 * eq_num + 1] * solInfoMap[n + 1].z2[i];
		}
		for( int i = 0; i < eq_num; ++i )
		{
			solInfoMap[n + 1].o[2 * eq_num + 2] += (*NtoOrt)( i ) * (*NtoOrt)( i );
		}
		solInfoMap[n + 1].o[2 * eq_num + 2] = sqrt( solInfoMap[n + 1].o[2 * eq_num + 2] );

		if( ( norm / solInfoMap[n + 1].o[2 * eq_num + 2] ).real() <= k11 )
		{
			for( int i = 0; i < eq_num; ++i )
			{
				(*NtoOrt)( i ) = (*NtoOrt)( i )	/ solInfoMap[n + 1].o[2 * eq_num + 2];
				solInfoMap[n + 1].z3[i] = (*NtoOrt)( i );
			}
		}
		else
		{
			for( int i = 0; i < eq_num; ++i )
			{
				omega2[ 2 * eq_num + 0 ] += (*NtoOrt)( i ) * solInfoMap[n + 1].z1[i];
			}
			for( int i = 0; i < eq_num; ++i )
			{
				(*NtoOrt)( i ) -= omega2[ 2 * eq_num + 0 ] * solInfoMap[n + 1].z1[i];
			}
			for( int i = 0; i < eq_num; ++i )
			{
				omega2[ 2 * eq_num + 1 ] += (*NtoOrt)( i ) * solInfoMap[n + 1].z2[i];
			}
			for( int i = 0; i < eq_num; ++i )
			{
				(*NtoOrt)( i ) -= omega2[ 2 * eq_num + 1 ] * solInfoMap[n + 1].z2[i];
			}
			for( int i = 0; i < eq_num; ++i )
			{
				omega2[ 2 * eq_num + 2 ] += (*NtoOrt)( i ) * (*NtoOrt)( i );
			}
			omega2[ 2 * eq_num + 2 ] = sqrt( omega2[ 2 * eq_num + 2 ] );

			for( int i = 0; i < eq_num; ++i )
			{
				solInfoMap[n + 1].z3[i] = (*NtoOrt)( i ) / omega2[ 2 * eq_num + 2 ];
				(*NtoOrt)( i ) = solInfoMap[n + 1].z3[i];
			}
			solInfoMap[n + 1].o[2 * eq_num + 0] += omega2[ 2 * eq_num + 0 ];
			solInfoMap[n + 1].o[2 * eq_num + 1] += omega2[ 2 * eq_num + 1 ];
			solInfoMap[n + 1].o[2 * eq_num + 2] = omega2[ 2 * eq_num + 2 ];
		}
	}
	else if( baseV == 4 )
	{
		for( int i = 0; i < eq_num; ++i )
		{
			norm += (*NtoOrt)( i ) * (*NtoOrt)( i );
			solInfoMap[n + 1].o[3 * eq_num + 0] += (*NtoOrt)( i ) * solInfoMap[n + 1].z1[i];
		}
		norm = sqrt( norm );
		for( int i = 0; i < eq_num; ++i )
		{
			(*NtoOrt)( i ) -= solInfoMap[n + 1].o[3 * eq_num + 0] * solInfoMap[n + 1].z1[i];
		}
		for( int i = 0; i < eq_num; ++i )
		{
			solInfoMap[n + 1].o[3 * eq_num + 1] += (*NtoOrt)( i ) * solInfoMap[n + 1].z2[i];
		}
		for( int i = 0; i < eq_num; ++i )
		{
			(*NtoOrt)( i ) -= solInfoMap[n + 1].o[3 * eq_num + 1] * solInfoMap[n + 1].z2[i];
		}
		for( int i = 0; i < eq_num; ++i )
		{
			solInfoMap[n + 1].o[3 * eq_num + 2] += (*NtoOrt)( i ) * solInfoMap[n + 1].z3[i];
		}
		for( int i = 0; i < eq_num; ++i )
		{
			(*NtoOrt)( i ) -= solInfoMap[n + 1].o[3 * eq_num + 2] * solInfoMap[n + 1].z3[i];
		}
		for( int i = 0; i < eq_num; ++i )
		{
			solInfoMap[n + 1].o[3 * eq_num + 3] += (*NtoOrt)( i ) * (*NtoOrt)( i );
		}
		solInfoMap[n + 1].o[3 * eq_num + 3] = sqrt( solInfoMap[n + 1].o[3 * eq_num + 3] );
		
		if( ( norm / solInfoMap[n + 1].o[3 * eq_num + 3] ).real() <= k11 )
		{
			for( int i = 0; i < eq_num; ++i )
			{
				(*NtoOrt)( i ) = (*NtoOrt)( i )	/ solInfoMap[n + 1].o[3 * eq_num + 3];
				solInfoMap[n + 1].z4[i] = (*NtoOrt)( i );
			}
		}
		else
		{
			//cout << "N4 > k11\n";
			for( int i = 0; i < eq_num; ++i )
			{
				omega2[ 3 * eq_num + 0 ] += (*NtoOrt)( i ) * solInfoMap[n + 1].z1[i];
			}
			for( int i = 0; i < eq_num; ++i )
			{
				(*NtoOrt)( i ) -= omega2[ 3 * eq_num + 0 ] * solInfoMap[n + 1].z1[i];
			}
			for( int i = 0; i < eq_num; ++i )
			{
				omega2[ 3 * eq_num + 1 ] += (*NtoOrt)( i ) * solInfoMap[n + 1].z2[i];
			}
			for( int i = 0; i < eq_num; ++i )
			{
				(*NtoOrt)( i ) -= omega2[ 3 * eq_num + 1 ] * solInfoMap[n + 1].z2[i];
			}
			for( int i = 0; i < eq_num; ++i )
			{
				omega2[ 3 * eq_num + 2 ] += (*NtoOrt)( i ) * solInfoMap[n + 1].z3[i];
			}
			for( int i = 0; i < eq_num; ++i )
			{
				(*NtoOrt)( i ) -= omega2[ 3 * eq_num + 2 ] * solInfoMap[n + 1].z3[i];
			}
			for( int i = 0; i < eq_num; ++i )
			{
				omega2[ 3 * eq_num + 3 ] += (*NtoOrt)( i ) * (*NtoOrt)( i );
			}
			omega2[ 3 * eq_num + 3 ] = sqrt( omega2[ 3 * eq_num + 3 ] );

			for( int i = 0; i < eq_num; ++i )
			{
				solInfoMap[n + 1].z4[i] = (*NtoOrt)( i ) / omega2[ 3 * eq_num + 3 ];
				(*NtoOrt)( i ) = solInfoMap[n + 1].z4[i];
			}
			solInfoMap[n + 1].o[3 * eq_num + 0] += omega2[ 3 * eq_num + 0 ];
			solInfoMap[n + 1].o[3 * eq_num + 1] += omega2[ 3 * eq_num + 1 ];
			solInfoMap[n + 1].o[3 * eq_num + 2] += omega2[ 3 * eq_num + 2 ];
			solInfoMap[n + 1].o[3 * eq_num + 3] = omega2[ 3 * eq_num + 3 ];
		}
	}
	else if( baseV == 5 )
	{
		for( int i = 0; i < eq_num; ++i )
		{
			solInfoMap[n + 1].o[4 * eq_num + 0] += (*NtoOrt)( i ) * solInfoMap[n + 1].z1[i];
			solInfoMap[n + 1].o[4 * eq_num + 1] += (*NtoOrt)( i ) * solInfoMap[n + 1].z2[i];
			solInfoMap[n + 1].o[4 * eq_num + 2] += (*NtoOrt)( i ) * solInfoMap[n + 1].z3[i];
			solInfoMap[n + 1].o[4 * eq_num + 3] += (*NtoOrt)( i ) * solInfoMap[n + 1].z4[i];
		}
		for( int i = 0; i < eq_num; ++i )
		{
			solInfoMap[n + 1].z5[i] = (*NtoOrt)( i ) - solInfoMap[n + 1].o[4 * eq_num + 0] * solInfoMap[n + 1].z1[i]
									- solInfoMap[n + 1].o[4 * eq_num + 1] * solInfoMap[n + 1].z2[i]
									- solInfoMap[n + 1].o[4 * eq_num + 2] * solInfoMap[n + 1].z3[i]
									- solInfoMap[n + 1].o[4 * eq_num + 3] * solInfoMap[n + 1].z4[i];
			(*NtoOrt)( i ) = solInfoMap[n + 1].z5[i];
		}
	}
}

template<class PL_NUM>
void OrthoBuilderGSh<PL_NUM>::buildSolution( vector<VarVect<PL_NUM> >* _mesh )
{
	vector<PL_NUM> M;
	vector<PL_NUM> f11;

	vector<PL_NUM> L;
	vector<PL_NUM> R;

	vector<PL_NUM> x1;

	int msize = eq_num / 2;				//caution here!
	M.resize( msize * msize, 0.0 );
	f11.resize( msize, 0.0 );

	L.resize( msize * msize, 0.0 );
	R.resize( msize * msize, 0.0 );

	x1.resize( msize, 0.0 );
	
	//simply supported plate NO CURRENT PASSING THROUGH THE BOUNDARY
	M[0 * msize + 0] = solInfoMap[Km - 1].z1[0];
	M[1 * msize + 0] = solInfoMap[Km - 1].z1[1];
	M[2 * msize + 0] = solInfoMap[Km - 1].z1[5];
	M[3 * msize + 0] = solInfoMap[Km - 1].z1[6];

	M[0 * msize + 1] = solInfoMap[Km - 1].z2[0];
	M[1 * msize + 1] = solInfoMap[Km - 1].z2[1];
	M[2 * msize + 1] = solInfoMap[Km - 1].z2[5];
	M[3 * msize + 1] = solInfoMap[Km - 1].z2[6];

	M[0 * msize + 2] = solInfoMap[Km - 1].z3[0];
	M[1 * msize + 2] = solInfoMap[Km - 1].z3[1];
	M[2 * msize + 2] = solInfoMap[Km - 1].z3[5];
	M[3 * msize + 2] = solInfoMap[Km - 1].z3[6];

	M[0 * msize + 3] = solInfoMap[Km - 1].z4[0];
	M[1 * msize + 3] = solInfoMap[Km - 1].z4[1];
	M[2 * msize + 3] = solInfoMap[Km - 1].z4[5];
	M[3 * msize + 3] = solInfoMap[Km - 1].z4[6];

	f11[0] = -solInfoMap[Km - 1].z5[0];
	f11[1] = -solInfoMap[Km - 1].z5[1];
	f11[2] = -solInfoMap[Km - 1].z5[5];
	f11[3] = -solInfoMap[Km - 1].z5[6];

	int posI = 0; 
	int posJ = 0;
	for( int i = 0; i < msize; ++i )
	{
		int found = 0;
		//for(int j = 0; j < msize; ++j )
		//{
			int j =0;
			if( fabs( M[i * msize + j].real() ) >= 0.0000000001 )
			{
				found = 1;
				posI = i;
				posJ = j;
				break;
			}
		//}
		if( found == 1 )
		{
			break;
		}
	}

	PL_NUM m0, m1, m2, m3;

	m0 = M[posI * msize + 0];
	m1 = M[posI * msize + 1];
	m2 = M[posI * msize + 2];
	m3 = M[posI * msize + 3];
	M[posI * msize + 0] = M[0 * msize + 0];
	M[posI * msize + 1] = M[0 * msize + 1];
	M[posI * msize + 2] = M[0 * msize + 2];
	M[posI * msize + 3] = M[0 * msize + 3];
	M[0 * msize + 0] = m0;
	M[0 * msize + 1] = m1;
	M[0 * msize + 2] = m2;
	M[0 * msize + 3] = m3;
	
	//m0 = M[0 * msize + posJ];
	//m1 = M[0 * msize + posJ];
	//m2 = M[0 * msize + posJ];
	//m3 = M[0 * msize + posJ];
	//M[0 * msize + posJ] = M[0 * msize + 0];
	//M[1 * msize + posJ] = M[1 * msize + 0];
	//M[2 * msize + posJ] = M[2 * msize + 0];
	//M[3 * msize + posJ] = M[3 * msize + 0];
	//M[0 * msize + 0] = m0;
	//M[1 * msize + 0] = m1;
	//M[2 * msize + 0] = m2;
	//M[3 * msize + 0] = m3;

	m0 = f11[posI];
	f11[posI] = f11[0];
	f11[0] = m0;

	for( int i = 0; i < msize; ++i )
	{
		L[i * msize + 0] = M[i * msize + 0];
		for( int j = i + 1; j < msize; ++j )		//TODO I think this is redundant
		{
			L[i * msize + j] = 0.0;
			R[j * msize + i] = 0.0;
		}
	}
	for( int i = 1; i < msize; ++i )
	{
		R[0 * msize + i] = M[0 * msize + i] / L[0 * msize + 0];
	}
	for( int i = 0; i < msize; ++i )
	{
		R[i * msize + i] = 1.0;
	}
	for( int i = 1; i < msize; ++i )
	{
		for( int k = i; k < msize; ++k )
		{
			L[k * msize + i] = M[k * msize + i];
			for( int l = 0; l < i; ++l )
			{
				L[k * msize + i] -= L[k * msize + l] * R[l * msize + i];
			}
		}
		for( int j = i + 1; j < msize; ++j )
		{
			R[i * msize + j] = M[i * msize + j];
			for( int l = 0; l < i; ++l )
			{
				R[i * msize + j] -= L[i * msize + l] * R[l * msize + j];
			}
			R[i * msize + j] /= L[i * msize + i];
		}
	}
	f11[0] = f11[0] / L[0 * msize + 0];
	for( int i = 1; i < msize; ++i )
	{
		for( int j = 0; j < i; ++j )
		{
			f11[i] -= L[i * msize + j] * f11[j];
		}
		f11[i] /= L[i * msize + i];
	}
	x1[msize - 1] = f11[msize - 1];
	for( int i = msize - 2; i >= 0; --i )
	{
		x1[i] = f11[i];
		for( int j = i + 1; j < msize; ++j )
		{
			x1[i] -= R[i * msize + j] * x1[j];
		}
	}

	for( int i = 0; i < msize; ++i )
	{
		solInfoMap[Km - 1].C[i] = x1[i];
	}

	//m0 = solInfoMap[Km - 1].C[0];								//we changed the order in M matrix, so we must change it in C too
	//solInfoMap[Km - 1].C[0] = solInfoMap[Km - 1].C[posJ];
	//solInfoMap[Km - 1].C[posJ] = m0;

	for( int j = 0; j < eq_num; ++j )
	{
		(*_mesh)[Km - 1].Nk1[j] = solInfoMap[Km - 1].C[0] * solInfoMap[Km - 1].z1[j]
						+ solInfoMap[Km - 1].C[1] * solInfoMap[Km - 1].z2[j]
						+ solInfoMap[Km - 1].C[2] * solInfoMap[Km - 1].z3[j]
						+ solInfoMap[Km - 1].C[3] * solInfoMap[Km - 1].z4[j]
						+ solInfoMap[Km - 1].z5[j];
	}

	for( int x = Km - 2; x >= 0; --x )			//calculate the coeeficients at all the other points and restore the solution there
	{
		solInfoMap[x].C[3] = ( solInfoMap[x + 1].C[3] 
								- solInfoMap[x + 1].o[4 * eq_num + 3] ) 
								/ solInfoMap[x + 1].o[3 * eq_num + 3];
		solInfoMap[x].C[2] = ( solInfoMap[x + 1].C[2] 
								- solInfoMap[x + 1].o[4 * eq_num + 2]
								- solInfoMap[x + 1].o[3 * eq_num + 2] * solInfoMap[x].C[3] ) 
								/ solInfoMap[x + 1].o[2 * eq_num + 2];
		solInfoMap[x].C[1] = ( solInfoMap[x + 1].C[1] 
								- solInfoMap[x + 1].o[4 * eq_num + 1]
								- solInfoMap[x + 1].o[3 * eq_num + 1] * solInfoMap[x].C[3]
								- solInfoMap[x + 1].o[2 * eq_num + 1] * solInfoMap[x].C[2] )
								/ solInfoMap[x + 1].o[1 * eq_num + 1];
		solInfoMap[x].C[0] = ( solInfoMap[x + 1].C[0] 
								- solInfoMap[x + 1].o[4 * eq_num + 0]
								- solInfoMap[x + 1].o[3 * eq_num + 0] * solInfoMap[x].C[3]
								- solInfoMap[x + 1].o[2 * eq_num + 0] * solInfoMap[x].C[2]
								- solInfoMap[x + 1].o[1 * eq_num + 0] * solInfoMap[x].C[1] )
								/ solInfoMap[x + 1].o[0 * eq_num + 0];
		for( int j = 0; j < eq_num; ++j )
		{
			(*_mesh)[x].Nk1[j] = solInfoMap[x].C[0] * solInfoMap[x].z1[j]
							+ solInfoMap[x].C[1] * solInfoMap[x].z2[j]
							+ solInfoMap[x].C[2] * solInfoMap[x].z3[j]
							+ solInfoMap[x].C[3] * solInfoMap[x].z4[j]
							+ solInfoMap[x].z5[j];
		}
	}
}

#endif