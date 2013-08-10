#ifndef _PLATE_1D_ADAMSBASHFORTH_
#define _PLATE_1D_ADAMSBASHFORTH_ 1

#include "plate_var_types.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <complex>

using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::complex;

class ABPrevPhis
{
public:
	ABPrevPhis();
	~ABPrevPhis();
	vector<PL_NUM> F_1;		//Function Phi of the system for basis vector of the -1 spatial step
	vector<PL_NUM> F_2;		//Function Phi of the system for basis vector of the -2 spatial step
	vector<PL_NUM> F_3;		//Function Phi of the system for basis vector of the -3 spatial step
	vector<PL_NUM> F;		//Function Phi of the system for basis vector of the current spatial step
};

class AdamsBashforth
{
public:
	AdamsBashforth();
	~AdamsBashforth();

	vector<ABPrevPhis> prevPhis;		//sets of 4 previous functions of the system for each of basis vectors 

	void calc( const vector<PL_NUM>& A, const vector<PL_NUM>& f, PL_NUM dx, int n, int hom, vector<PL_NUM>* x );	//hom == 0 if homogeneous
	void setEqNum( int _eq_num );
private:
	int eq_num;
};

#endif