#include "VarVect.h"

VarVectAdj::VarVectAdj()
{

}

VarVectAdj::~VarVectAdj()
{

}

VarVectAdj::VarVectAdj( int _eq_num )
{
	N.resize( _eq_num, 0.0 );

	N1.resize( _eq_num, 0.0 );
	d1N1.resize( _eq_num, 0.0 );
	d2N1.resize( _eq_num, 0.0 );
}

void VarVectAdj::setup( int _eq_num )
{
	N.resize( _eq_num, 0.0 );

	N1.resize( _eq_num, 0.0 );
	d1N1.resize( _eq_num, 0.0 );
	d2N1.resize( _eq_num, 0.0 );
}
