#ifndef _PLATE_1D_HYPERDUAL_
#define _PLATE_1D_HYPERDUAL_ 1

#include <vector>
#include <iostream>

using std::vector;
using std::ostream;

template<class D_PRES, int NN>
class HPD;

//--------------------------------------------------------------------
// addition

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator+( const HPD<D_PRES, NN>& lhs, const HPD<D_PRES, NN>& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i] + rhs.elems[i];
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator+( const D_PRES& lhs, const HPD<D_PRES, NN>& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = rhs.elems[i];
	}
	ret.elems[0] += lhs;
	return ret;
}

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator+( const HPD<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i];
	}
	ret.elems[0] += rhs;
	return ret;
}

//----------------------------------------------------------------
// subtraction

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator-( const HPD<D_PRES, NN>& lhs, const HPD<D_PRES, NN>& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i] - rhs.elems[i];
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator-( const D_PRES& lhs, const HPD<D_PRES, NN>& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = -rhs.elems[i];
	}
	ret.elems[0] += lhs;
	return ret;
}

template<class D_PRES, int NN>
const HPD<D_PRES, NN> operator-( const HPD<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	HPD<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i];
	}
	ret.elems[0] -= rhs;
	return ret;
}

//----------------------------------------------------------------
// cout

template<class D_PRES, int NN>
ostream& operator<<( ostream& os, const HPD<D_PRES, NN>& item )
{
	for( int i = 0; i < NN; ++i )
	{
		os << item.elems[i] << " ";
	}
	os << item.elems[NN];
	return os;
}

//----------------------------------------------------------------
// class itself

template<class D_PRES, int NN>
class HPD		//HPD = hyper dual
{
public:
	vector<D_PRES> elems;

	HPD() {	elems.resize( NN + 1, 0.0 ); }
	HPD( D_PRES a ) {	elems.resize( NN + 1, 0.0 ); elems[0] = a; }
	D_PRES real() { return elems[0]; }

	template<class D_PRES_, int NN_>
	friend ostream& operator<<( ostream& os, const HPD& item );
	
	HPD& operator=( const HPD& rhs );
	HPD& operator=( const D_PRES& rhs );

	friend const HPD operator+<>( const HPD& lhs, const HPD& rhs );
	friend const HPD operator+<>( const D_PRES& lhs, const HPD& rhs );
	friend const HPD operator+<>( const HPD& lhs, const D_PRES& rhs );

	friend const HPD operator-<>( const HPD& lhs, const HPD& rhs );
	friend const HPD operator-<>( const D_PRES& lhs, const HPD& rhs );
	friend const HPD operator-<>( const HPD& lhs, const D_PRES& rhs );

private:
	HPD( const HPD& rhs );
};

template<class D_PRES, int NN>
HPD<D_PRES, NN>::HPD( const HPD& rhs )
{
	cout << " d-----\n";
}

template<class D_PRES, int NN>
HPD<D_PRES, NN>& HPD<D_PRES, NN>::operator=( const HPD& rhs )
{
	if( this != &rhs )
	{
		for( int i = 0; i < NN + 1; ++i )
		{
			elems[i] = rhs.elems[i];
		}
	}
	return *this;
}

template<class D_PRES, int NN>
HPD<D_PRES, NN>& HPD<D_PRES, NN>::operator=( const D_PRES& rhs )
{
	for( int i = 0; i < NN + 1; ++i )
	{
		elems[i] = 0.0;
	}
	elems[0] = rhs;
	return *this;
}

#endif