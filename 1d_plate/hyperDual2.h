#ifndef _PLATE_1D_HYPERDUAL2_
#define _PLATE_1D_HYPERDUAL2_ 1

#include <vector>
#include <iostream>

using std::vector;
using std::ostream;

template<class D_PRES, int NN>
class HPD2;

//--------------------------------------------------------------------
// addition

template<class D_PRES, int NN>
const HPD2<D_PRES, NN> operator+( const HPD2<D_PRES, NN>& lhs, const HPD2<D_PRES, NN>& rhs )
{
	HPD2<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i] + rhs.elems[i];
	}
	for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
	{
		ret.elems2[i] = lhs.elems2[i] + rhs.elems2[i];
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD2<D_PRES, NN> operator+( const D_PRES& lhs, const HPD2<D_PRES, NN>& rhs )
{
	HPD2<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = rhs.elems[i];
	}
	for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
	{
		ret.elems2[i] = rhs.elems2[i];
	}
	ret.elems[0] += lhs;
	return ret;
}

template<class D_PRES, int NN>
const HPD2<D_PRES, NN> operator+( const HPD2<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	HPD2<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i];
	}
	for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
	{
		ret.elems2[i] = lhs.elems2[i];
	}
	ret.elems[0] += rhs;
	return ret;
}

//----------------------------------------------------------------
// subtraction

template<class D_PRES, int NN>
const HPD2<D_PRES, NN> operator-( const HPD2<D_PRES, NN>& lhs, const HPD2<D_PRES, NN>& rhs )
{
	HPD2<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i] - rhs.elems[i];
	}
	for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
	{
		ret.elems2[i] = lhs.elems2[i] - rhs.elems2[i];
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD2<D_PRES, NN> operator-( const D_PRES& lhs, const HPD2<D_PRES, NN>& rhs )
{
	HPD2<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = -rhs.elems[i];
	}
	for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
	{
		ret.elems2[i] = -rhs.elems2[i];
	}
	ret.elems[0] += lhs;
	return ret;
}

template<class D_PRES, int NN>
const HPD2<D_PRES, NN> operator-( const HPD2<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	HPD2<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i];
	}
	for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
	{
		ret.elems2[i] = lhs.elems2[i];
	}
	ret.elems[0] -= rhs;
	return ret;
}

template<class D_PRES, int NN>
const HPD2<D_PRES, NN> operator-( const HPD2<D_PRES, NN>& rhs )
{
	HPD2<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = -rhs.elems[i];
	}
	for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
	{
		ret.elems2[i] = -rhs.elems2[i];
	}
	return ret;
}

//--------------------------------------------------------------------
// multiplication

template<class D_PRES, int NN>
const HPD2<D_PRES, NN> operator*( const HPD2<D_PRES, NN>& lhs, const HPD2<D_PRES, NN>& rhs )
{
	HPD2<D_PRES, NN> ret;

	ret.elems[0] = lhs.elems[0] * rhs.elems[0];
	for( int i = 1; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[0] * rhs.elems[i] + lhs.elems[i] * rhs.elems[0];
	}
	int ind = 0;
	for( int i = 1; i < NN + 1; ++i )
	{
		for( int j = i; j < NN + 1; ++j )
		{
			ind = ( i - 1 ) * NN - ( i - 2 ) * ( i - 1 ) / 2 + ( j - i );
			ret.elems2[ind] = lhs.elems[0] * rhs.elems2[ind] + lhs.elems[i] * rhs.elems[j] + lhs.elems[j] * rhs.elems[i] + lhs.elems2[ind] * rhs.elems[0];
		}
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD2<D_PRES, NN> operator*( const D_PRES& lhs, const HPD2<D_PRES, NN>& rhs )
{
	HPD2<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = rhs.elems[i] * lhs;
	}
	for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
	{
		ret.elems2[i] = rhs.elems2[i] * lhs;
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD2<D_PRES, NN> operator*( const HPD2<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	HPD2<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i] * rhs;
	}
	for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
	{
		ret.elems2[i] = lhs.elems2[i] * rhs;
	}
	return ret;
}

//--------------------------------------------------------------------
// division

template<class D_PRES, int NN>
const HPD2<D_PRES, NN> operator/( const HPD2<D_PRES, NN>& lhs, const HPD2<D_PRES, NN>& rhs )
{
	HPD2<D_PRES, NN> ret;

	ret.elems[0] = lhs.elems[0] / rhs.elems[0];
	for( int i = 1; i < NN + 1; ++i )
	{
		ret.elems[i] = ( lhs.elems[i] * rhs.elems[0] - lhs.elems[0] * rhs.elems[i] ) / ( rhs.elems[0] * rhs.elems[0] ); 
	}
	int ind = 0;
	for( int i = 1; i < NN + 1; ++i )
	{
		for( int j = i; j < NN + 1; ++j )
		{
			ind = ( i - 1 ) * NN - ( i - 2 ) * ( i - 1 ) / 2 + ( j - i );
			ret.elems2[ind] = lhs.elems[0] * ( 2.0l / rhs.elems[0] / rhs.elems[0] / rhs.elems[0] * rhs.elems[i] * rhs.elems[j] - 1.0l / rhs.elems[0] / rhs.elems[0] * rhs.elems2[ind] ) 
							- lhs.elems[i] / rhs.elems[0] / rhs.elems[0] * rhs.elems[j]
							- lhs.elems[j] / rhs.elems[0] / rhs.elems[0] * rhs.elems[i]
							+ lhs.elems2[ind] / rhs.elems[0];
		}
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD2<D_PRES, NN> operator/( const D_PRES& lhs, const HPD2<D_PRES, NN>& rhs )
{
	HPD2<D_PRES, NN> ret;

	ret.elems[0] = lhs / rhs.elems[0];
	for( int i = 1; i < NN + 1; ++i )
	{
		ret.elems[i] = ( -lhs * rhs.elems[i] ) / ( rhs.elems[0] * rhs.elems[0] ); 
	}
	int ind = 0;
	for( int i = 1; i < NN + 1; ++i )
	{
		for( int j = i; j < NN + 1; ++j )
		{
			ind = ( i - 1 ) * NN - ( i - 2 ) * ( i - 1 ) / 2 + ( j - i );
			ret.elems2[ind] = lhs * ( 2.0l / rhs.elems[0] / rhs.elems[0] / rhs.elems[0] * rhs.elems[i] * rhs.elems[j] - 1.0l / rhs.elems[0] / rhs.elems[0] * rhs.elems2[ind] );
		}
	}
	return ret;
}

template<class D_PRES, int NN>
const HPD2<D_PRES, NN> operator/( const HPD2<D_PRES, NN>& lhs, const D_PRES& rhs )
{
	HPD2<D_PRES, NN> ret;

	for( int i = 0; i < NN + 1; ++i )
	{
		ret.elems[i] = lhs.elems[i] / rhs;
	}
	for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
	{
		ret.elems2[i] = lhs.elems2[i] / rhs;
	}
	return ret;
}

//----------------------------------------------------------------
// cout

template<class D_PRES, int NN>
ostream& operator<<( ostream& os, const HPD2<D_PRES, NN>& item )
{
	for( int i = 0; i < NN; ++i )
	{
		os << item.elems[i] << " ";
	}
	os << item.elems[NN] << endl;
	for( int i = 0; i < NN; ++i )
	{
		os << item.elems2[i];
		int ind = 0;
		for( int j = 1; j <= i; ++j )
		{
			ind += NN - j;
			os << " " << item.elems2[i + ind];
		}
		os << endl;
	}
	return os;
}

//----------------------------------------------------------------
// class itself

template<class D_PRES, int NN>
class HPD2		//HPD = hyper dual
{
public:
	//vector<D_PRES> elems;
	D_PRES elems[NN + 1];	//contains the actual value and first order partial derivatives
	D_PRES elems2[( NN + 1 ) * NN / 2];		//containd the second order derivatives - the matrix is symmetric => we need only "half" of it

	HPD2() 
	{	
		//for( int i = 0; i < NN + 1; ++i )
		//{
		//	elems[i] = 0.0;
		//	elems2[i] = 0.0;
		//}
		//for( int i = NN + 1; i < ( NN + 1 ) * NN / 2; ++i )
		//{
		//	elems2[i] = 0.0;
		//}
	}
	//	elems.resize( NN + 1, 0.0 ); }
	HPD2( D_PRES a )
	//{	elems.resize( NN + 1, 0.0 ); elems[0] = a; }
	{
		for( int i = 0; i < NN + 1; ++i )
		{
			elems[i] = 0.0;
			elems2[i] = 0.0;
		}
		for( int i = NN + 1; i < ( NN + 1 ) * NN / 2; ++i )
		{
			elems2[i] = 0.0;
		}
		elems[0] = a;
	}
	D_PRES real() const { return elems[0]; }

	template<class D_PRES_, int NN_>
	friend ostream& operator<<( ostream& os, const HPD2& item );
	
	HPD2& operator=( const HPD2& rhs );
	HPD2& operator=( const D_PRES& rhs );

	HPD2& operator+=( const HPD2& rhs );
	HPD2& operator+=( const D_PRES& rhs );

	HPD2& operator-=( const HPD2& rhs );
	HPD2& operator-=( const D_PRES& rhs );

	HPD2& operator*=( const HPD2& rhs );
	HPD2& operator*=( const D_PRES& rhs );

	HPD2& operator/=( const HPD2& rhs );
	HPD2& operator/=( const D_PRES& rhs );

	friend const HPD2 operator+<>( const HPD2& lhs, const HPD2& rhs );
	friend const HPD2 operator+<>( const D_PRES& lhs, const HPD2& rhs );
	friend const HPD2 operator+<>( const HPD2& lhs, const D_PRES& rhs );

	friend const HPD2 operator-<>( const HPD2& lhs, const HPD2& rhs );
	friend const HPD2 operator-<>( const D_PRES& lhs, const HPD2& rhs );
	friend const HPD2 operator-<>( const HPD2& lhs, const D_PRES& rhs );

	friend const HPD2 operator-<>( const HPD2& rhs );		//negate

	friend const HPD2 operator*<>( const HPD2& lhs, const HPD2& rhs );
	friend const HPD2 operator*<>( const D_PRES& lhs, const HPD2& rhs );
	friend const HPD2 operator*<>( const HPD2& lhs, const D_PRES& rhs );

	friend const HPD2 operator/<>( const HPD2& lhs, const HPD2& rhs );
	friend const HPD2 operator/<>( const D_PRES& lhs, const HPD2& rhs );
	friend const HPD2 operator/<>( const HPD2& lhs, const D_PRES& rhs );

//private:
	//HPD( const HPD2& rhs );
};

//template<class D_PRES, int NN>
//HPD2<D_PRES, NN>::HPD( const HPD2& rhs )
//{
//	cout << " d-----\n";
//}

template<class D_PRES, int NN>
HPD2<D_PRES, NN>& HPD2<D_PRES, NN>::operator=( const HPD2& rhs )
{
	if( this != &rhs )
	{
		for( int i = 0; i < NN + 1; ++i )
		{
			elems[i] = rhs.elems[i];
		}
		for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
		{
			elems2[i] = rhs.elems2[i];
		}
	}
	return *this;
}

template<class D_PRES, int NN>
HPD2<D_PRES, NN>& HPD2<D_PRES, NN>::operator=( const D_PRES& rhs )
{
	for( int i = 0; i < NN + 1; ++i )
	{
		elems[i] = 0.0;
	}
	for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
	{
		elems2[i] = 0.0;
	}
	elems[0] = rhs;
	return *this;
}

template<class D_PRES, int NN>
HPD2<D_PRES, NN>& HPD2<D_PRES, NN>::operator+=( const HPD2& rhs )
{
	if( this != &rhs )
	{
		for( int i = 0; i < NN + 1; ++i )
		{
			elems[i] += rhs.elems[i];
		}
		for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
		{
			elems2[i] += rhs.elems2[i];
		}
	}
	return *this;
}

template<class D_PRES, int NN>
HPD2<D_PRES, NN>& HPD2<D_PRES, NN>::operator+=( const D_PRES& rhs )
{
	elems[0] += rhs;
	return *this;
}

template<class D_PRES, int NN>
HPD2<D_PRES, NN>& HPD2<D_PRES, NN>::operator-=( const HPD2& rhs )
{
	if( this != &rhs )
	{
		for( int i = 0; i < NN + 1; ++i )
		{
			elems[i] -= rhs.elems[i];
		}
		for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
		{
			elems2[i] -= rhs.elems2[i];
		}
	}
	return *this;
}

template<class D_PRES, int NN>
HPD2<D_PRES, NN>& HPD2<D_PRES, NN>::operator-=( const D_PRES& rhs )
{
	elems[0] -= rhs;
	return *this;
}

template<class D_PRES, int NN>
HPD2<D_PRES, NN>& HPD2<D_PRES, NN>::operator*=( const HPD2& rhs )
{
	//HPD2<D_PRES, NN> tmp;
	//tmp = ( *this ) * rhs;
	( *this ) = ( *this ) * rhs;

	return *this;
}

template<class D_PRES, int NN>
HPD2<D_PRES, NN>& HPD2<D_PRES, NN>::operator*=( const D_PRES& rhs )
{
	for( int i = 0; i < NN; ++i )
	{
		elems[i] *= rhs;
	}
	for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
	{
		elems2[i] *= rhs;
	}
	return *this;
}

template<class D_PRES, int NN>
HPD2<D_PRES, NN>& HPD2<D_PRES, NN>::operator/=( const HPD2& rhs )
{
	( *this ) = ( *this ) / rhs;

	return *this;
}

template<class D_PRES, int NN>
HPD2<D_PRES, NN>& HPD2<D_PRES, NN>::operator/=( const D_PRES& rhs )
{
	for( int i = 0; i < NN; ++i )
	{
		elems[i] /= rhs;
	}
	for( int i = 0; i < ( NN + 1 ) * NN / 2; ++i )
	{
		elems2[i] /= rhs;
	}
	return *this;
}

template<class D_PRES, int NN>
HPD2<D_PRES, NN> sin( const HPD2<D_PRES, NN>& arg )
{
	HPD2<D_PRES, NN> ret;
	D_PRES x0 = arg.elems[0];
	ret = sin( x0 ) + cos( x0 ) * ( arg - x0 );
	int ind = 0;
	for( int i = 1; i < NN + 1; ++i )
	{
		for( int j = i; j < NN + 1; ++j )
		{
			ind = ( i - 1 ) * NN - ( i - 2 ) * ( i - 1 ) / 2 + ( j - i );
			ret.elems2[ind] = -sin( x0 ) * arg.elems[i] * arg.elems[j] + cos( x0 ) * arg.elems2[ind];
		}
	}
	return ret;
}

template<class D_PRES, int NN>
HPD2<D_PRES, NN> cos( const HPD2<D_PRES, NN>& arg )
{
	HPD2<D_PRES, NN> ret;
	D_PRES x0 = arg.elems[0];
	ret = cos( x0 ) - sin( x0 ) * ( arg - x0 );
	int ind = 0;
	for( int i = 1; i < NN + 1; ++i )
	{
		for( int j = i; j < NN + 1; ++j )
		{
			ind = ( i - 1 ) * NN - ( i - 2 ) * ( i - 1 ) / 2 + ( j - i );
			ret.elems2[ind] = -cos( x0 ) * arg.elems[i] * arg.elems[j] - sin( x0 ) * arg.elems2[ind];
		}
	}
	return ret;
}

template<class D_PRES, int NN>
HPD2<D_PRES, NN> exp( const HPD2<D_PRES, NN>& arg )
{
	HPD2<D_PRES, NN> ret;
	D_PRES x0 = arg.elems[0];
	ret = exp( x0 ) + exp( x0 ) * ( arg - x0 ) + exp( x0 ) * ( arg - x0 ) * ( arg - x0 ) / 2.0l;
	return ret;
}

template<class D_PRES, int NN>
HPD2<D_PRES, NN> sqrt( const HPD2<D_PRES, NN>& arg )
{
	HPD2<D_PRES, NN> ret;
	D_PRES x0 = arg.elems[0];
	ret = sqrt( x0 ) + 0.5l / sqrt( x0 ) * ( arg - x0 ) - 1.0l / 8.0l / sqrt( x0 ) / sqrt( x0 ) / sqrt( x0 ) * ( arg - x0 ) * ( arg - x0 );
	return ret;
}

#endif