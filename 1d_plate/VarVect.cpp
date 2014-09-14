#include "VarVect.h"

VarVectAdj::VarVectAdj()
{

}

VarVectAdj::~VarVectAdj()
{

}

VarVectAdj::VarVectAdj( int _eq_num )
{
	Nk.resize( _eq_num, 0.0 );

	Nk0.resize( _eq_num, 0.0 );
	d1N0.resize( _eq_num, 0.0 );
	d2N0.resize( _eq_num, 0.0 );
}

void VarVectAdj::setup( int _eq_num )
{
	Nk.resize( _eq_num, 0.0 );

	Nk0.resize( _eq_num, 0.0 );
	d1N0.resize( _eq_num, 0.0 );
	d2N0.resize( _eq_num, 0.0 );
}
