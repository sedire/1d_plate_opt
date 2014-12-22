#ifndef _PLATE_1D_PARAMSACK_
#define _PLATE_1D_PARAMSACK_ 1

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include "plate_var_types.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;
using std::stringstream;

class ParamSack	//this one is used to store all the parameters for the program
{
public:
	ParamSack();
	~ParamSack();

	void loadFromFile( string inputFname );

private:
	N_PRES E1;				//Young's modulus
	N_PRES E2;				//Young's modulus
	N_PRES nu21;			//Poisson's ratio	
	N_PRES rho;				//composite's density
	N_PRES sigma_x;			//electric conductivity
	N_PRES mu;
	N_PRES eps_0;
	N_PRES eps_x;

	N_PRES h;				//thickness of the plate
	N_PRES a;				//width of the plate

	int currentType;
	int numberOfCurrentParams;
	vector<N_PRES> currentParams;

	int stressType;
	int numberOfStressParams;
	vector<N_PRES> stressParams;

	N_PRES totalTimeT1;
	N_PRES switchTimeT0;

	//parameters if the numerical methods
	int nodeNumOnY;
	N_PRES dt;
	int maxNewtonIter;
	N_PRES orthonormCheckEps;
	N_PRES newtonEps;
	N_PRES newtonAlmostZero;
};

#endif