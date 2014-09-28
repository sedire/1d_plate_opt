#ifndef _PLATE_1D_ADJSOLVER_
#define _PLATE_1D_ADJSOLVER_ 1

#include "plate_var_types.h" 
#include <iostream>
#include "SolverPar.h"
#include "VarVect.h"
#include <Eigen\Eigen>
#include "RungeKutta.h"
#include "OrthoBuilder.h"
#include <sstream>
#include <fstream>

using std::cout;
using std::endl;
using std::stringstream;
using std::ofstream;
using namespace Eigen;

class AdjSolver		
{
public:
	AdjSolver();
	~AdjSolver();
	void loadParamsFromStruct( const SolverPar& loadFrom );
	void setPrimalSolnData( N_PRES* _primSoln );
	void setAdjointSolnData( N_PRES* _adjSoln );

	N_PRES doStep();		//the adjoint problem is solved backwards in time
	void decreaseTime();
	N_PRES getCurTime();
	int getCurTimeStep();
	void dumpSol( int fNum );

private:
//parameters of the material
	N_PRES E1;				//Young's modulus
	N_PRES E2;				//Young's modulus
	N_PRES nu21;			//Poisson's ratio	
	N_PRES rho;				//composite's density

	N_PRES sigma_x;			//electric conductivity
	N_PRES sigma_x_mu;

	N_PRES h;				//thickness of the plate
	N_PRES a;				//width of the plate
//scheme params
	N_PRES dt;
	const int totTimeSteps;
	int Km;
	N_PRES dx;
//time track
	N_PRES curTime;
	int curTimeStep;
	N_PRES switchTime;
//stress
	N_PRES J0;
	N_PRES J0_1;
	N_PRES tauSin;
	N_PRES tauSin_1;
	N_PRES tauExp;
	N_PRES tauExp_1;

	N_PRES p0;				//constant mechanical load
	N_PRES tauP;
	N_PRES rad;

	int stress_type;
	int current_type;
//Newmark params
	N_PRES beta;
//some other
	N_PRES B11;
	N_PRES B22;
	N_PRES B12;
	N_PRES By0;
	N_PRES By1;                                      // in considered boundary-value problem

	N_PRES eps_0;
	N_PRES eps_x;
	N_PRES eps_x_0;

//other
	RungeKutta<N_PRES>* rungeKutta;
	OrthoBuilder<N_PRES>* orthoBuilder;

	const int eq_num;
	vector<VarVectAdj> mesh;
	vector<N_PRES> newmarkA;
	vector<N_PRES> newmarkB;

	Matrix<N_PRES, EQ_NUM, EQ_NUM> matrA;
	Matrix<N_PRES, EQ_NUM, 1> vectF;

	N_PRES* primSoln;	//pointer to where the primal solution is stored
	N_PRES* primSolnDt;	//pointer to the time derivative approximation of the primal solution -- computed here
	N_PRES* adjSoln;	//here we will store the adjoint solution

//basis vectors for superposition part
	Matrix<N_PRES, EQ_NUM, 1> N1;
	Matrix<N_PRES, EQ_NUM, 1> N2;
	Matrix<N_PRES, EQ_NUM, 1> N3;
	Matrix<N_PRES, EQ_NUM, 1> N4;
	Matrix<N_PRES, EQ_NUM, 1> N5;

	void calcNewmarkAB( int y );
	void calcSystemMatrices( int y );
};

#endif