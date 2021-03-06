#ifndef _PLATE_1D_VARVECT_
#define _PLATE_1D_VARVECT_ 1

#include <vector>
#include "plate_var_types.h"
#include <complex>

using std::vector;
using std::complex;

class VarVect
{
public:
	VarVect();
	VarVect( int _eq_num );
	~VarVect();
	void setup( int _eq_num );

	vector<PL_NUM> Nk;
	vector<PL_NUM> Nk1;
	vector<PL_NUM> d1N;
	vector<PL_NUM> d2N;

	vector<PL_NUM> Nk0;			//don't really know why we need these. some computational tricks, probably
	vector<PL_NUM> d1N0;
	vector<PL_NUM> d2N0;
};

#endif