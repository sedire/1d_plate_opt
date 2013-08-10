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

class OrthoBuilder
{
public:
	OrthoBuilder();
	virtual ~OrthoBuilder();
	virtual void setParams( int _Km );
	//virtual void orthonorm( int baseV, int n, vector<PL_NUM>* NtoOrt ) {};
	virtual void orthonorm( int baseV, int n, Matrix<PL_NUM, EQ_NUM, 1>* NtoOrt ) {};
	virtual void buildSolution( vector<VarVect>* _mesh ) {};
	virtual void flushO( int x );
	virtual void setInitVects( const Matrix<PL_NUM, EQ_NUM, 1> &N1, const Matrix<PL_NUM, EQ_NUM, 1> &N2, const Matrix<PL_NUM, EQ_NUM, 1> &N3, const Matrix<PL_NUM, EQ_NUM, 1> &N4, const Matrix<PL_NUM, EQ_NUM, 1> &N5 );
	virtual void orthogTest( int x );
protected:
	int eq_num;
	int Km;
	vector<SolInfo> solInfoMap;
};

class OrthoBuilderGodunov : public OrthoBuilder
{
public:
	OrthoBuilderGodunov() {};
	~OrthoBuilderGodunov() {};
	//void orthonorm( int baseV, int n, vector<PL_NUM>* NtoOrt );
	//void orthonorm( int baseV, int n, PL_NUM* NtoOrt );
	void buildSolution( vector<VarVect>* _mesh );
};

class OrthoBuilderGSh : public OrthoBuilder
{
public:
	OrthoBuilderGSh() {};
	~OrthoBuilderGSh() {};
	//void orthonorm( int baseV, int n, vector<PL_NUM>* NtoOrt );
	void orthonorm( int baseV, int n, Matrix<PL_NUM, EQ_NUM, 1>* NtoOrt );
	void buildSolution( vector<VarVect>* _mesh );
};

#endif