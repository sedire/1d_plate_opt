#ifndef _PLATE_1D_ORTHOBUILDER_
#define _PLATE_1D_ORTHOBUILDER_ 1

#include "plate_var_types.h"
#include "VarVect.h"
#include <iostream>
#include <vector>
#include <Eigen/Eigen>
#include <complex>
#include "time.h"

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

	//vector<PL_NUM> C;

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
	inline virtual void orthonorm( int baseV, int n, Matrix<PL_NUM, EQ_NUM, 1>* NtoOrt ) {};
	inline virtual int checkOrtho( int n, const Matrix<PL_NUM, EQ_NUM, 1>& N2orthog, 
											const Matrix<PL_NUM, EQ_NUM, 1>& N3orthog,
											const Matrix<PL_NUM, EQ_NUM, 1>& N4orthog, 
											const Matrix<PL_NUM, EQ_NUM, 1>& N5orthog, 

											const Matrix<PL_NUM, EQ_NUM, 1>& N2orig,
											const Matrix<PL_NUM, EQ_NUM, 1>& N3orig,
											const Matrix<PL_NUM, EQ_NUM, 1>& N4orig,
											const Matrix<PL_NUM, EQ_NUM, 1>& N5orig ) { return 1; };
	inline virtual void setNextSolVects( int n, const Matrix<PL_NUM, EQ_NUM, 1>& N1, 
									const Matrix<PL_NUM, EQ_NUM, 1>& N2, 
									const Matrix<PL_NUM, EQ_NUM, 1>& N3,
									const Matrix<PL_NUM, EQ_NUM, 1>& N4, 
									const Matrix<PL_NUM, EQ_NUM, 1>& N5 ) {};
	inline virtual void normNextSolVects( int n, Matrix<PL_NUM, EQ_NUM, 1>* N1, 
									Matrix<PL_NUM, EQ_NUM, 1>* N2, 
									Matrix<PL_NUM, EQ_NUM, 1>* N3,
									Matrix<PL_NUM, EQ_NUM, 1>* N4 ) {};
	virtual void buildSolution( vector<VarVect<PL_NUM> >* _mesh ) {};
	virtual void buildSolutionAdj( vector<VarVectAdj>* _mesh ) {};
	virtual void flushO( int x );
	virtual void setInitVects( const Matrix<PL_NUM, EQ_NUM, 1> &N1, const Matrix<PL_NUM, EQ_NUM, 1> &N2, const Matrix<PL_NUM, EQ_NUM, 1> &N3, const Matrix<PL_NUM, EQ_NUM, 1> &N4, const Matrix<PL_NUM, EQ_NUM, 1> &N5 );
	virtual void orthogTest( int x );
	virtual void resetOrthoDoneInfo();
	virtual void setOrthoDoneInfo( int y );

	time_t orthoStart;
	time_t orthoTotal;
protected:
	int eq_num;
	int Km;
	vector<SolInfo<PL_NUM> > solInfoMap;
	vector<PL_NUM> omega2;
	vector<bool> orthoDone;
	vector<PL_NUM> Cx1;		//to keep coefficients for the superposition solution from the next orthonormalization interval
	vector<PL_NUM> Cx;		//to keep coefficients for the superposition solution from the current orthonormalization interval
};

template<class PL_NUM>
class OrthoBuilderGSh : public OrthoBuilder<PL_NUM>
{
public:
	OrthoBuilderGSh() {};
	~OrthoBuilderGSh() {};
	//void orthonorm( int baseV, int n, vector<PL_NUM>* NtoOrt );
	inline void orthonorm( int baseV, int n, Matrix<PL_NUM, EQ_NUM, 1>* NtoOrt );
	inline int checkOrtho( int n, const Matrix<PL_NUM, EQ_NUM, 1>& N2, 
											const Matrix<PL_NUM, EQ_NUM, 1>& N3,
											const Matrix<PL_NUM, EQ_NUM, 1>& N4, 
											const Matrix<PL_NUM, EQ_NUM, 1>& N5, 

											const Matrix<PL_NUM, EQ_NUM, 1>& N2orig,
											const Matrix<PL_NUM, EQ_NUM, 1>& N3orig,
											const Matrix<PL_NUM, EQ_NUM, 1>& N4orig,
											const Matrix<PL_NUM, EQ_NUM, 1>& N5orig );
	inline void setNextSolVects( int n, const Matrix<PL_NUM, EQ_NUM, 1>& N1, 
									const Matrix<PL_NUM, EQ_NUM, 1>& N2, 
									const Matrix<PL_NUM, EQ_NUM, 1>& N3,
									const Matrix<PL_NUM, EQ_NUM, 1>& N4, 
									const Matrix<PL_NUM, EQ_NUM, 1>& N5 );
	inline void normNextSolVects( int n, Matrix<PL_NUM, EQ_NUM, 1>* N1, 
									Matrix<PL_NUM, EQ_NUM, 1>* N2, 
									Matrix<PL_NUM, EQ_NUM, 1>* N3,
									Matrix<PL_NUM, EQ_NUM, 1>* N4 );
	void buildSolution( vector<VarVect<PL_NUM> >* _mesh );
	void buildSolutionAdj( vector<VarVectAdj>* _mesh );
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

	//C.resize( EQ_NUM / 2 );
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
OrthoBuilder<PL_NUM>::OrthoBuilder() :
	eq_num( EQ_NUM ),
	orthoStart( 0 ),
	orthoTotal( 0 )
{
	omega2.resize( eq_num * eq_num, 0.0 );
	Cx1.resize( eq_num / 2, 0.0 );
	Cx.resize( eq_num / 2, 0.0 );
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

	orthoDone.resize( Km - 1, 0 );
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
void OrthoBuilder<PL_NUM>::resetOrthoDoneInfo()
{
	for( int y = 0; y < orthoDone.size(); ++y )
	{
		orthoDone[y] = false;
	}
}

template<class PL_NUM>
void OrthoBuilder<PL_NUM>::setOrthoDoneInfo( int y )
{
	orthoDone[y] = true;
}

template<class PL_NUM>
int OrthoBuilderGSh<PL_NUM>::checkOrtho( int n, const Matrix<PL_NUM, EQ_NUM, 1>& N2orthog, 
											const Matrix<PL_NUM, EQ_NUM, 1>& N3orthog,
											const Matrix<PL_NUM, EQ_NUM, 1>& N4orthog, 
											const Matrix<PL_NUM, EQ_NUM, 1>& N5orthog, 

											const Matrix<PL_NUM, EQ_NUM, 1>& N2orig,
											const Matrix<PL_NUM, EQ_NUM, 1>& N3orig,
											const Matrix<PL_NUM, EQ_NUM, 1>& N4orig,
											const Matrix<PL_NUM, EQ_NUM, 1>& N5orig )
{
	int ret = 0;
	N_PRES eps = ORTHONORM_CHECK_EPS;
	
	if( N2orthog.lpNorm<Infinity>() * solInfoMap[n + 1].o[1 * eq_num + 1] < eps * N2orig.lpNorm<Infinity>() ||
		N3orthog.lpNorm<Infinity>() * solInfoMap[n + 1].o[2 * eq_num + 2] < eps * N3orig.lpNorm<Infinity>() || 
		N4orthog.lpNorm<Infinity>() * solInfoMap[n + 1].o[3 * eq_num + 3] < eps * N4orig.lpNorm<Infinity>() ||
		N5orthog.lpNorm<Infinity>() < eps * N5orig.lpNorm<Infinity>() )
	{
		ret = 1;
		orthoDone[n] = true;
	}

	return ret;
}

template<class PL_NUM>
inline void OrthoBuilderGSh<PL_NUM>::setNextSolVects( int n, const Matrix<PL_NUM, EQ_NUM, 1>& N1, 
									const Matrix<PL_NUM, EQ_NUM, 1>& N2, 
									const Matrix<PL_NUM, EQ_NUM, 1>& N3,
									const Matrix<PL_NUM, EQ_NUM, 1>& N4, 
									const Matrix<PL_NUM, EQ_NUM, 1>& N5 )
{
	for( int i = 0; i < eq_num; ++i )
	{
		solInfoMap[n + 1].z1[i] = N1( i );
		solInfoMap[n + 1].z2[i] = N2( i );
		solInfoMap[n + 1].z3[i] = N3( i );
		solInfoMap[n + 1].z4[i] = N4( i );
		solInfoMap[n + 1].z5[i] = N5( i );
	}
}

template<class PL_NUM>
inline void OrthoBuilderGSh<PL_NUM>::normNextSolVects( int n, Matrix<PL_NUM, EQ_NUM, 1>* N1, 
									Matrix<PL_NUM, EQ_NUM, 1>* N2, 
									Matrix<PL_NUM, EQ_NUM, 1>* N3,
									Matrix<PL_NUM, EQ_NUM, 1>* N4 )
{
	(*N1) /= solInfoMap[n + 1].o[0 * eq_num + 0];
	(*N2) /= solInfoMap[n + 1].o[1 * eq_num + 1];
	(*N3) /= solInfoMap[n + 1].o[2 * eq_num + 2];
	(*N4) /= solInfoMap[n + 1].o[3 * eq_num + 3];
}

template<class PL_NUM>
void OrthoBuilderGSh<PL_NUM>::orthonorm( int baseV, int n, Matrix<PL_NUM, EQ_NUM, 1>* NtoOrt )
{
	orthoStart = time( 0 );

	long double k11 = 1.414213562373095;
	PL_NUM norm = 0;
	for( int i = 0; i < eq_num * eq_num; ++i )
	{
		omega2[i] = 0.0;
	}

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

		if( norm / solInfoMap[n + 1].o[1 * eq_num + 1] <= k11 )
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

		if( norm / solInfoMap[n + 1].o[2 * eq_num + 2] <= k11 )
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
		
		if( norm / solInfoMap[n + 1].o[3 * eq_num + 3] <= k11 )
		{
			for( int i = 0; i < eq_num; ++i )
			{
				(*NtoOrt)( i ) = (*NtoOrt)( i )	/ solInfoMap[n + 1].o[3 * eq_num + 3];
				solInfoMap[n + 1].z4[i] = (*NtoOrt)( i );
			}
		}
		else
		{
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
			(*NtoOrt)( i ) = (*NtoOrt)( i ) - solInfoMap[n + 1].o[4 * eq_num + 0] * solInfoMap[n + 1].z1[i]
									- solInfoMap[n + 1].o[4 * eq_num + 1] * solInfoMap[n + 1].z2[i]
									- solInfoMap[n + 1].o[4 * eq_num + 2] * solInfoMap[n + 1].z3[i]
									- solInfoMap[n + 1].o[4 * eq_num + 3] * solInfoMap[n + 1].z4[i];
			solInfoMap[n + 1].z5[i] = (*NtoOrt)( i );
		}
	}

	orthoTotal += time( 0 ) - orthoStart;
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
			if( fabs( M[i * msize + j] ) >= 0.0000000001l )
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
		Cx1[i] = x1[i];
	}

	//m0 = solInfoMap[Km - 1].C[0];								//we changed the order in M matrix, so we must change it in C too
	//solInfoMap[Km - 1].C[0] = solInfoMap[Km - 1].C[posJ];
	//solInfoMap[Km - 1].C[posJ] = m0;

	for( int j = 0; j < eq_num; ++j )
	{
		(*_mesh)[Km - 1].Nk1[j] = Cx1[0] * solInfoMap[Km - 1].z1[j]
						+ Cx1[1] * solInfoMap[Km - 1].z2[j]
						+ Cx1[2] * solInfoMap[Km - 1].z3[j]
						+ Cx1[3] * solInfoMap[Km - 1].z4[j]
						+ solInfoMap[Km - 1].z5[j];

	}

	for( int x = Km - 2; x >= 0; --x )			//calculate the coeeficients at all the other points and restore the solution there
	{
		if( orthoDone[x] == true )
		{
			Cx[3] = ( Cx1[3] 
						- solInfoMap[x + 1].o[4 * eq_num + 3] ) 
						/ solInfoMap[x + 1].o[3 * eq_num + 3];
			Cx[2] = ( Cx1[2] 
						- solInfoMap[x + 1].o[4 * eq_num + 2]
						- solInfoMap[x + 1].o[3 * eq_num + 2] * Cx[3] ) 
						/ solInfoMap[x + 1].o[2 * eq_num + 2];
			Cx[1] = ( Cx1[1] 
						- solInfoMap[x + 1].o[4 * eq_num + 1]
						- solInfoMap[x + 1].o[3 * eq_num + 1] * Cx[3]
						- solInfoMap[x + 1].o[2 * eq_num + 1] * Cx[2] )
						/ solInfoMap[x + 1].o[1 * eq_num + 1];
			Cx[0] = ( Cx1[0] 
						- solInfoMap[x + 1].o[4 * eq_num + 0]
						- solInfoMap[x + 1].o[3 * eq_num + 0] * Cx[3]
						- solInfoMap[x + 1].o[2 * eq_num + 0] * Cx[2]
						- solInfoMap[x + 1].o[1 * eq_num + 0] * Cx[1] )
						/ solInfoMap[x + 1].o[0 * eq_num + 0];

			//update the coefficients for the next interval
			for( int i = 0; i < eq_num / 2; ++i )
			{
				Cx1[i] = Cx[i];
			}
		}
		else
		{
			for( int i = 0; i < eq_num / 2; ++i )
			{
				Cx[i] = Cx1[i];
			}
			//... and there is no need to update the coefficients for the next interval
		}
		for( int j = 0; j < eq_num; ++j )
		{
			(*_mesh)[x].Nk1[j] = Cx[0] * solInfoMap[x].z1[j]
							+ Cx[1] * solInfoMap[x].z2[j]
							+ Cx[2] * solInfoMap[x].z3[j]
							+ Cx[3] * solInfoMap[x].z4[j]
							+ solInfoMap[x].z5[j];
		}
	}
}

template<class PL_NUM>
void OrthoBuilderGSh<PL_NUM>::buildSolutionAdj( vector<VarVectAdj>* _mesh )
{
}

template<>
void inline OrthoBuilderGSh<N_PRES>::buildSolutionAdj( vector<VarVectAdj>* _mesh )
{
	static const int msize = EQ_NUM / 2;
	Matrix<N_PRES, msize, msize, RowMajor> M;
	Matrix<N_PRES, msize, 1> f11;
	Matrix<N_PRES, msize, 1> x1;
	Matrix<N_PRES, msize, 1> res;
	Matrix<N_PRES, msize, 1> dx;

	for( int i = 0; i < msize; ++i )
	{
		for( int j = 0; j < msize; ++j )
		{
			M( j, i ) = 0.0;
		}
		f11( i ) = 0.0;
		x1( i ) = 0.0;
	}
	
	//simply supported plate NO CURRENT PASSING THROUGH THE BOUNDARY
	M( 0, 0 ) = solInfoMap[Km - 1].z1[2];
	M( 1, 0 ) = solInfoMap[Km - 1].z1[3];
	M( 2, 0 ) = solInfoMap[Km - 1].z1[4];
	M( 3, 0 ) = solInfoMap[Km - 1].z1[7];

	M( 0, 1 ) = solInfoMap[Km - 1].z2[2];
	M( 1, 1 ) = solInfoMap[Km - 1].z2[3];
	M( 2, 1 ) = solInfoMap[Km - 1].z2[4];
	M( 3, 1 ) = solInfoMap[Km - 1].z2[7];

	M( 0, 2 ) = solInfoMap[Km - 1].z3[2];
	M( 1, 2 ) = solInfoMap[Km - 1].z3[3];
	M( 2, 2 ) = solInfoMap[Km - 1].z3[4];
	M( 3, 2 ) = solInfoMap[Km - 1].z3[7];

	M( 0, 3 ) = solInfoMap[Km - 1].z4[2];
	M( 1, 3 ) = solInfoMap[Km - 1].z4[3];
	M( 2, 3 ) = solInfoMap[Km - 1].z4[4];
	M( 3, 3 ) = solInfoMap[Km - 1].z4[7];

	f11( 0 ) = -solInfoMap[Km - 1].z5[2];
	f11( 1 ) = -solInfoMap[Km - 1].z5[3];
	f11( 2 ) = -solInfoMap[Km - 1].z5[4];
	f11( 3 ) = -solInfoMap[Km - 1].z5[7];

	x1 = M.fullPivLu().solve( f11 );
	res = f11 - M * x1;
	//cout << " --------- " << res.lpNorm<Infinity>() / x1.lpNorm<Infinity>() << endl;

	for( int i = 0; i < msize; ++i )
	{
		Cx1[i] = x1( i );
	}

	//m0 = solInfoMap[Km - 1].C[0];								//we changed the order in M matrix, so we must change it in C too
	//solInfoMap[Km - 1].C[0] = solInfoMap[Km - 1].C[posJ];
	//solInfoMap[Km - 1].C[posJ] = m0;

	for( int j = 0; j < eq_num; ++j )
	{
		(*_mesh)[Km - 1].N[j] = Cx1[0] * solInfoMap[Km - 1].z1[j]
						+ Cx1[1] * solInfoMap[Km - 1].z2[j]
						+ Cx1[2] * solInfoMap[Km - 1].z3[j]
						+ Cx1[3] * solInfoMap[Km - 1].z4[j]
						+ solInfoMap[Km - 1].z5[j];
	}

	for( int x = Km - 2; x >= 0; --x )			//calculate the coeeficients at all the other points and restore the solution there
	{
		if( orthoDone[x] == true )
		{
			Cx[3] = ( Cx1[3] 
						- solInfoMap[x + 1].o[4 * eq_num + 3] ) 
						/ solInfoMap[x + 1].o[3 * eq_num + 3];
			Cx[2] = ( Cx1[2] 
						- solInfoMap[x + 1].o[4 * eq_num + 2]
						- solInfoMap[x + 1].o[3 * eq_num + 2] * Cx[3] ) 
						/ solInfoMap[x + 1].o[2 * eq_num + 2];
			Cx[1] = ( Cx1[1] 
						- solInfoMap[x + 1].o[4 * eq_num + 1]
						- solInfoMap[x + 1].o[3 * eq_num + 1] * Cx[3]
						- solInfoMap[x + 1].o[2 * eq_num + 1] * Cx[2] )
						/ solInfoMap[x + 1].o[1 * eq_num + 1];
			Cx[0] = ( Cx1[0] 
						- solInfoMap[x + 1].o[4 * eq_num + 0]
						- solInfoMap[x + 1].o[3 * eq_num + 0] * Cx[3]
						- solInfoMap[x + 1].o[2 * eq_num + 0] * Cx[2]
						- solInfoMap[x + 1].o[1 * eq_num + 0] * Cx[1] )
						/ solInfoMap[x + 1].o[0 * eq_num + 0];

			//update the coefficients for the next interval
			for( int i = 0; i < eq_num / 2; ++i )
			{
				Cx1[i] = Cx[i];
			}
		}
		else
		{
			for( int i = 0; i < eq_num / 2; ++i )
			{
				Cx[i] = Cx1[i];
			}
			//... and there is no need to update the coefficients for the next interval
		}
		for( int j = 0; j < eq_num; ++j )
		{
			(*_mesh)[x].N[j] = Cx[0] * solInfoMap[x].z1[j]
							+ Cx[1] * solInfoMap[x].z2[j]
							+ Cx[2] * solInfoMap[x].z3[j]
							+ Cx[3] * solInfoMap[x].z4[j]
							+ solInfoMap[x].z5[j];
		}
	}
}

#endif