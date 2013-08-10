#include "Optimizer.h"

Optimizer::Optimizer( Solver* _solver, N_PRES _weight, N_PRES _J0h, N_PRES _tauh, N_PRES _charTime ):
	solver( _solver ),
	weight( _weight ),
	J0h( _J0h ),
	tauh( _tauh ),
	charTime( _charTime )
{
	parVect( 0 ) = 0.0;
	parVect( 1 ) = 0.0;

	parVectPrev( 0 ) = 0.0;
	parVectPrev( 1 ) = 0.0;

	objVal = 0.0;
	objValPrev = 0.0;

	gk1( 0 ) = 1.0;
	gk1( 1 ) = 1.0;

	gk( 0 ) = 1.0;
	gk( 1 ) = 1.0;

	dk1( 0 ) = 0.0;
	dk1( 1 ) = 0.0;

	dk( 0 ) = 0.0;
	dk( 1 ) = 0.0;

	alphak = 0.0;

	wolfeDelta = 0.0001;
	wolfeSigma = 0.99;

	epsK = 0.000001;
	thetaUpdate = 0.5;

	phiPrAA = 0.0;
	phiPrBB = 0.0;

	gammaLineS = 0.66;

	dmpo.open( "optSteps.txt" );
	dmpo.close();
	dmpo.open( "optSteps.txt", ofstream::app );
}

void Optimizer::optimize( N_PRES Jstart, N_PRES tauStart )
{
	int count = 0;
	const N_PRES threshold = 0.000000000000001;
	parVect( 0 ) = Jstart;
	parVect( 1 ) = tauStart;

	while( gk1.norm() >= threshold )
	{
		cout << " optimization step " << count << endl;

		calc1stOrdOptInfo( parVect( 0 ), parVect( 1 ), &objVal, &gk1 );
		if( count != 0 )
		{
			dk1 = -gk1 + calcBettaN_() * dk;
		}
		else
		{
			dk1 = -gk1;
		}

		dmpo << parVect( 0 ) << " " << parVect( 1 ) << " " << objVal << " " << gk1( 0 ) << " " << gk1( 1 ) << " " << dk1( 0 ) << " " << dk1( 1 ) << " " << calcBettaN_() << " " << gk1.norm() << " " << threshold << endl;
		cout << " objVal " << objVal << endl;
		cout << " grad " << gk1 << endl;
		cout << " dir " << dk1 << endl;
		cout << " betta " << calcBettaN_() << endl;

		if( gk1.norm() < threshold )
		{
			break;
		}

		phi0 = objVal;
		phiPr0 = gk1.transpose() * dk1;
		if( lineSearch() != 1 )
		{
			cout << " optimization failed in line search\n";
			return;
		}

		parVectPrev = parVect;
		objValPrev = objVal;
		dk = dk1;
		gk = gk1;

		parVect = parVect  + linSb * dk1; 

		++count;
	}
}

int Optimizer::lineSearch()
{
	if( gk1.transpose() * dk1 < 0 )
	{
		phiA = phi0;
		phiPrA = gk1.transpose() * dk1;
	}
	else
	{
		cout << " something weird happened: (gk1, dk1) >= 0\n";
		return 0;
	}
	Matrix<N_PRES, 2, 1> newPar;
	N_PRES al = 0.1;
	newPar = parVect + al * dk1;

	while( newPar( 1 ) > 0.0064 || newPar( 1 ) < 0.0048 || 
		newPar( 0 ) <= 0.0 )
	{
		al /= 2.0;
		newPar = parVect + al * dk1;
	}
	
	N_PRES objValFora = 0.0;
	Matrix<N_PRES, 2, 1> gradFora;

	N_PRES objValForb = 0.0;
	Matrix<N_PRES, 2, 1> gradForb;
	al /= 1.1;
	do
	{
		al *= 1.1;
		newPar = parVect + al * dk1;
		calc1stOrdOptInfo( newPar( 0 ), newPar( 1 ), &objValForb, &gradForb );

		cout << " ===\n";
	}while( gradForb.transpose() * dk1 < 0 );

	if( gradForb.transpose() * dk1 < 0 )
	{
		cout << " something weird happened: (gradForb, dk1) < 0\n";
		return 0;
	}

	linSa = 0.0;
	linSb = al;
	phiB = objValForb;
	phiPrB = gradForb.transpose() * dk1;

	//initial a, phi(a), b, phi(b) are found satisfying phi(a) <= phi(0) + eps_k, phi'(a) < 0, phi'(b) >= 0. As in Hager's paper.

	while( !( ( phiB - phi0 <= wolfeDelta * linSb * gk1.transpose() * dk1 && 
		gradForb.transpose() * dk1 >= wolfeSigma * gk1.transpose() * dk1 ) ||
		( ( 2.0 * wolfeDelta - 1 ) * phiPr0 >= phiPrB &&
		phiPrB >= wolfeSigma * phiPr0 &&
		phiB <= phi0 + epsK ) ) )		//repeat until we satisfy either wolfe conditions of approximate wolfe conditions
	{
		cout << " ~~~\n";
		N_PRES newa = 0.0;
		N_PRES newb = 0.0;
		secant2( &newa, & newb );
		if( newb - newa > gammaLineS * ( linSb - linSa ) )
		{
			N_PRES newc = ( newa + newb ) / 2.0;
			N_PRES newaa = 0.0;
			N_PRES newbb = 0.0;
			update( newa, newb, newc, &newaa, &newbb );
			newa = newaa;
			newb = newbb;
		}

		linSa = newa;
		linSb = newb;

		newPar = parVect + linSb * dk1;
		calc1stOrdOptInfo( newPar( 0 ), newPar( 1 ), &objValForb, &gradForb );
		phiB = objValForb;
		phiPrB = gradForb.transpose() * dk1;

		newPar = parVect + linSa * dk1;
		calc1stOrdOptInfo( newPar( 0 ), newPar( 1 ), &objValFora, &gradFora );
		phiA = objValFora;
		phiPrA = gradFora.transpose() * dk1;
	}
	return 1;
}

void Optimizer::secant2( N_PRES* _a, N_PRES* _b )
{
	linSc = ( linSa * phiPrB - linSb * phiPrA ) / ( phiPrB - phiPrA );		//TODO check it
	N_PRES A = 0.0;
	N_PRES B = 0.0;
	update( linSa, linSb, linSc, &A, &B );

	N_PRES _c = 0.0;
	if( linSc == B )
	{
		_c = ( linSb * phiPrC - linSc * phiPrB ) / ( phiPrC - phiPrB );
	}
	if( linSc == A )
	{
		_c = ( linSa * phiPrC - linSc * phiPrA ) / ( phiPrC - phiPrA );
	}
	if( linSc == A || linSc == B )
	{
		update( A, B, _c, _a, _b );
	}
	else
	{
		*_a = A;
		*_b = B;
	}
	return;
}

void Optimizer::update( N_PRES a, N_PRES b, N_PRES c, N_PRES* _a, N_PRES* _b )
{
	if( c <= a || c >= b )
	{
		*_a = a;
		*_b = b;
		return;
	}

	N_PRES objValForc;
	Matrix<N_PRES, 2, 1> gradForc;
	Matrix<N_PRES, 2, 1> newParc = parVect + c * dk1;

	calc1stOrdOptInfo( newParc( 0 ), newParc( 1 ), &objValForc, &gradForc );
	phiC = objValForc;
	phiPrC = gradForc.transpose() * dk1;

	if( phiPrC >= 0.0 )
	{
		*_a = a;
		*_b = c;
		return;
	}
	if( phiPrC < 0.0 && phiC <= phi0 + epsK )
	{
		*_a = c;
		*_b = b;

		return;
	}
	if( phiPrC < 0.0 && phiC > phi0 + epsK )
	{
		N_PRES aa = a;
		N_PRES bb = c;

		while( 1 )
		{
			linSd = ( 1.0 - thetaUpdate ) * aa + thetaUpdate * bb;

			N_PRES objValFord;
			Matrix<N_PRES, 2, 1> gradFord;
			Matrix<N_PRES, 2, 1> newPard = parVect + linSd * dk1;
			calc1stOrdOptInfo( newPard( 0 ), newPard( 1 ), &objValFord, &gradFord );
			phiD = objValFord;
			phiPrD = gradFord.transpose() * dk1;
			if( phiPrD >= 0.0 )
			{
				*_b = linSd;
				*_a = aa;
				break;
			}
			if( phiPrD < 0.0 && phiD <= phi0 + epsK )
			{
				aa = linSd;
				continue;
			}
			if( phiPrD < 0.0 && phiD > phi0 + epsK )
			{
				bb = linSd;
				continue;
			}
		}
		return;
	}
}

void Optimizer::calc1stOrdOptInfo( long double J0, long double tau, long double* _objVal, Matrix<N_PRES, 2, 1>* _gk )
{
	PL_NUM J0begin;
	PL_NUM tauBegin;

	for( int i = 0; i < 2; ++i )
	{
		if( i == 0 )
		{
			J0begin = PL_NUM( J0, J0h );
			tauBegin = PL_NUM( tau, 0.0 );
		}
		else if( i == 1 )
		{
			J0begin = PL_NUM( J0, 0.0 );
			tauBegin = PL_NUM( tau, tauh );
		}
		solver->setTask( J0begin, tauBegin );
		solver->calcConsts();
		PL_NUM funcVal = calcFuncVal();

		if( i == 0 )
		{
			( *_gk )( 0 ) = funcVal.imag() / J0h * J0_SCALE + weight;
		}
		else if( i == 1 )
		{
			( *_gk )( 1 ) = funcVal.imag() / tauh;
		}
		*_objVal = funcVal.real() + J0 * weight;
	}
}

PL_NUM Optimizer::calcFuncVal()
{
	PL_NUM funcVal( 0.0, 0.0 );

	while( solver->cur_t.real() <= charTime )
	{
		PL_NUM val;
		val = solver->do_step();
		funcVal += val * val;

		solver->cur_t += solver->dt;
		++( solver->curTimeStep );
	}

	return funcVal;
}

N_PRES Optimizer::calcBettaN()
{
	Matrix<N_PRES, 2, 1> yk = gk1 - gk;

	return 1.0 / ( dk.transpose() * yk ) * ( yk - 2.0 * dk * ( yk.transpose() * yk ) / ( dk.transpose() * yk ) ).transpose() * gk1;
}

N_PRES Optimizer::calcBettaN_()
{
	N_PRES oldBetta = calcBettaN();
	N_PRES etta = 0.01;		//as in Hager's paper
	N_PRES ettak = -1.0 / ( dk.norm() * min( etta, gk.norm() ) );
	
	return max( oldBetta, ettak );
}

N_PRES Optimizer::min( N_PRES a, N_PRES b )
{
	N_PRES ret = a;
	if( b < ret )
	{
		ret = b;
	}
	return ret;
}

N_PRES Optimizer::max( N_PRES a, N_PRES b )
{
	N_PRES ret = a;
	if( b > ret )
	{
		ret = b;
	}
	return ret;
}