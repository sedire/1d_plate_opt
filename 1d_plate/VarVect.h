#ifndef _PLATE_1D_VARVECT_
#define _PLATE_1D_VARVECT_ 1

#include <vector>
#include "plate_var_types.h"
#include <complex>

using std::vector;
using std::complex;

template<class PL_NUM>
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

	vector<PL_NUM> Nk0;
	vector<PL_NUM> d1N0;
	vector<PL_NUM> d2N0;
};


template<class PL_NUM>
VarVect<PL_NUM>::VarVect()
{

}

template<class PL_NUM>
VarVect<PL_NUM>::~VarVect()
{

}

template<class PL_NUM>
VarVect<PL_NUM>::VarVect( int _eq_num )
{
	Nk.resize( _eq_num, 0.0 );
	Nk1.resize( _eq_num, 0.0 );
	d1N.resize( _eq_num, 0.0 );
	d2N.resize( _eq_num, 0.0 );

	Nk0.resize( _eq_num, 0.0 );
	d1N0.resize( _eq_num, 0.0 );
	d2N0.resize( _eq_num, 0.0 );
}

template<class PL_NUM>
void VarVect<PL_NUM>::setup( int _eq_num )
{
	Nk.resize( _eq_num, 0.0 );
	Nk1.resize( _eq_num, 0.0 );
	d1N.resize( _eq_num, 0.0 );
	d2N.resize( _eq_num, 0.0 );

	Nk0.resize( _eq_num, 0.0 );
	d1N0.resize( _eq_num, 0.0 );
	d2N0.resize( _eq_num, 0.0 );
}

#endif