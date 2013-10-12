#ifndef _PLATE_1D_OPTIMIZER_
#define _PLATE_1D_OPTIMIZER_ 1

#include <iostream>
#include <fstream>
#include "Eigen/Eigen"
#include "plate_var_types.h"
#include "hyperDual.h"
#include "Solver.h"
#include "asa_user.h"
#include "HagerOptFuncs.h"

using namespace Eigen;
using std::cout;
using std::endl;
using std::ofstream;

//double calc1stOrdOptInfoCG_DES( double* g, double* x, long n );

template<class PL_NUM>
class Optimizer
{
public:
	Optimizer( Solver<PL_NUM>* _solver, N_PRES _weightJ, N_PRES _weightB , N_PRES charTime );
	void optimize( const Matrix<N_PRES, GRAD_SIZE, 1>& params );
	void optimizeCG_DES( const Matrix<N_PRES, GRAD_SIZE, 1>& params );
	void optimizeASA( const Matrix<N_PRES, GRAD_SIZE, 1>& params );
	void optimizeNewton( const Matrix<N_PRES, GRAD_SIZE, 1>& params );

private:
	void calc1stOrdOptInfo( const Matrix<N_PRES, GRAD_SIZE, 1>& curVal, long double* _objVal, Matrix<N_PRES, GRAD_SIZE, 1>* _gk );
	void calc2ndOrdOptInfo( const Matrix<N_PRES, GRAD_SIZE, 1>& curVal, long double* _objVal, Matrix<N_PRES, GRAD_SIZE, 1>* _gk, Matrix<N_PRES, GRAD_SIZE, GRAD_SIZE>* _hess );

	PL_NUM calcFuncVal();
	N_PRES calcBettaN();
	N_PRES calcBettaN_();
	int lineSearch();
	void secant2( N_PRES* _a, N_PRES* _b );
	void update( N_PRES a, N_PRES b, N_PRES c, N_PRES* _a, N_PRES* _b );

	N_PRES min( N_PRES a, N_PRES b );
	N_PRES max( N_PRES a, N_PRES b );

	Solver<PL_NUM>* solver;

	N_PRES weightJ;
	N_PRES weightB;

	N_PRES charTime;

	Matrix<N_PRES, GRAD_SIZE, 1> parVect;
	Matrix<N_PRES, GRAD_SIZE, 1> parVectPrev;

	N_PRES objVal;
	N_PRES objValPrev;
	
	Matrix<N_PRES, GRAD_SIZE, 1> gk1;		//new gradient
	Matrix<N_PRES, GRAD_SIZE, 1> gk;		//gradient from prev step

	Matrix<N_PRES, GRAD_SIZE, 1> dk1;		//new descent direction
	Matrix<N_PRES, GRAD_SIZE, 1> dk;		//descent direction from previous step

	N_PRES alphak;

	N_PRES linSa;
	N_PRES linSb;
	N_PRES linSc;
	N_PRES linSd;

	N_PRES phi0;
	N_PRES phiPr0;
	N_PRES phiAk;
	N_PRES phiA;
	N_PRES phiPrA;
	N_PRES phiB;
	N_PRES phiPrB;
	N_PRES phiC;
	N_PRES phiPrC;
	N_PRES phiD;
	N_PRES phiPrD;

	N_PRES phiPrAA;
	N_PRES phiPrBB;

	N_PRES wolfeDelta;
	N_PRES wolfeSigma;

	N_PRES epsK;
	N_PRES thetaUpdate;
	N_PRES gammaLineS;

	ofstream dmpo;
};

template<class PL_NUM>
Optimizer<PL_NUM>::Optimizer( Solver<PL_NUM>* _solver, N_PRES _weightJ, N_PRES _weightB, N_PRES _charTime ):
	solver( _solver ),
	weightJ( _weightJ ),
	weightB( _weightB ),
	charTime( _charTime )
{
	for( int i = 0; i < GRAD_SIZE; ++i )
	{
		parVect( i ) = 0.0;
		parVectPrev( i ) = 0.0;

		gk1( i ) = 1.0;
		gk( i ) = 1.0;
		dk1( i ) = 0.0;
		dk( i ) = 0.0;
	}

	objVal = 0.0;
	objValPrev = 0.0;

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

template<class PL_NUM>
void Optimizer<PL_NUM>::optimize( const Matrix<N_PRES, GRAD_SIZE, 1>& params )
{
	int count = 0;
	const N_PRES threshold = 0.000001;
	for( int i = 0; i < GRAD_SIZE; ++i )
	{
		parVect( i ) = params( i );
	}

	while( gk1.lpNorm<Infinity>() > threshold )
	{
		cout << " =======\n";
		cout << " optimization step " << count << endl;
		cout << " =======\n";

		calc1stOrdOptInfo( parVect, &objVal, &gk1 );

		//TODO: delete me
		//return;

		if( count != 0 )
		{
			dk1 = -gk1 + calcBettaN_() * dk;
		}
		else
		{
			dk1 = -gk1;
		}

		dmpo << " -- par: " << parVect( 0 ) << " " << parVect( 1 ) << " " << parVect( 2 ) << " -- obj: " << objVal; 
		dmpo << " -- j weight: " << parVect( 0 ) * parVect( 0 ) / 4.0l / ( 1.0l + M_PI * M_PI )
			* parVect( 1 ) * exp( -2.0l * charTime / parVect( 1 ) ) * ( M_PI * M_PI * ( exp( 2.0l * charTime / parVect( 1 ) ) - 1 )
			- M_PI * sin( 2.0l * M_PI * charTime / parVect( 1 ) ) + cos( 2.0l * M_PI * charTime / parVect( 1 ) ) - 1.0l );
		dmpo << " -- grad: "<< gk1( 0 ) << " " << gk1( 1 ) << " " << gk1( 2 ) << " -- dir: ";
		dmpo << dk1( 0 ) << " " << dk1( 1 ) << " " << dk1( 2 ) << " --- " << calcBettaN_() << " " << gk1.lpNorm<Infinity>() << " " << threshold << endl;

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

template<class PL_NUM>
void Optimizer<PL_NUM>::optimizeNewton( const Matrix<N_PRES, GRAD_SIZE, 1>& params )
{
	int count = 0;
	const N_PRES threshold = 0.000000000000001;

	Matrix<N_PRES, GRAD_SIZE, 1> gradient;
	for( int i = 0; i < GRAD_SIZE; ++i )
	{
		gradient( i ) = 1.0l;
		parVect( i ) = params( i );
	}

	Matrix<N_PRES, GRAD_SIZE, GRAD_SIZE> hessian;

	while( gradient.norm() >= threshold )
	{
		cout << " =======\n";
		cout << " optimization step " << count << endl;
		cout << " =======\n";

		calc2ndOrdOptInfo( parVect, &objVal, &gradient, &hessian );

		dmpo << " -- par: " << parVect( 0 ) << " " << parVect( 1 ) << " " << parVect( 2 ) << " -- obj: " << objVal << " -- grad: "; 
		dmpo << gradient( 0 ) << " " << gradient( 1 ) << " " << gradient( 2 ) << " -- dir: ";
		dmpo << " " << gradient.norm() << " " << threshold << endl;

		parVect -= 0.0001 * hessian.inverse() * gradient;

		++count;
	}
}

//template<class PL_NUM>
//void Optimizer<PL_NUM>::optimizeCG_DES( const Matrix<N_PRES, GRAD_SIZE, 1>& params )
//{
//	cout << "optimizeCG_DES enter\n";
//
//	const N_PRES threshold = 1.e-6;
//	double* x = new double[GRAD_SIZE];
//	for( int i = 0; i < GRAD_SIZE; ++i )
//	{
//		x[i] = params( i );
//	}
//
//	cg_parameter Parm;
//	cg_default( &Parm );    /* set default parameter values */
//    Parm.PrintLevel = 3;
//	Parm.PrintParms = 0;
//
//	cg_descent( x, GRAD_SIZE, 0, &Parm, threshold, calcValCG_DES, calcGradCG_DES, calc1stOrdOptInfoCG_DES, 0 );
//
//	delete[] x;
//}

template<class PL_NUM>
void Optimizer<PL_NUM>::optimizeASA( const Matrix<N_PRES, GRAD_SIZE, 1>& params )
{
	cout << "optimizeASA enter\n";

	const N_PRES threshold = 1.e-6;
	double* x = new double[GRAD_SIZE];
	for( int i = 0; i < GRAD_SIZE; ++i )
	{
		x[i] = params( i );
	}
	double* lo = new double[GRAD_SIZE];
	double* hi = new double[GRAD_SIZE];

	lo[0] = -1.0;
	lo[1] = 0.00001;
	lo[2] = 0.0;
	hi[0] = 1.0;
	hi[1] = 100.0;
	hi[2] = 1.0;

	asacg_parm cgParm;
    asa_parm asaParm;
	asa_cg_default( &cgParm );
    asa_default( &asaParm );
    cgParm.PrintParms = TRUE;
    cgParm.PrintLevel = 3;
    asaParm.PrintParms = TRUE;
    asaParm.PrintLevel = 3;

	//cg_descent( x, GRAD_SIZE, 0, &Parm, threshold, calcValCG_DES, calcGradCG_DES, calc1stOrdOptInfoCG_DES, 0 );
	asa_cg( x, lo, hi, GRAD_SIZE, NULL, &cgParm, &asaParm, threshold, calcValASA, calcGradASA, calc1stOrdOptInfoASA, 0, 0 ) ;

	cout << "\n\n===============\nASA optimization complete. X is:\n";
	for( int i = 0; i < GRAD_SIZE; ++i )
	{
		cout << x[i] << endl;
	}

	delete[] x;
	delete[] lo;
	delete[] hi;
}

template<class PL_NUM>
int Optimizer<PL_NUM>::lineSearch()
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
	Matrix<N_PRES, GRAD_SIZE, 1> newPar;
	N_PRES al = 0.1;
	N_PRES alStep = 1.4;
	newPar = parVect + al * dk1;

	//while( newPar( 1 ) > 0.0064 || newPar( 1 ) < 0.0048 || 
	//	newPar( 0 ) <= 0.0000000001 ||
	//	newPar( 2 ) <= 0.0000000001 || newPar( 2 ) >= 5.0 )
	//{
	//	al /= 2.0;
	//	newPar = parVect + al * dk1;
	//	cout << "lin search " << al << " " << newPar << endl;
	//}
	
	N_PRES objValFora = 0.0;
	Matrix<N_PRES, GRAD_SIZE, 1> gradFora;

	N_PRES objValForb = 0.0;
	Matrix<N_PRES, GRAD_SIZE, 1> gradForb;

	//N_PRES oldObjValForb = 0.0;
	//Matrix<N_PRES, GRAD_SIZE, 1> oldGradForb;

	int it = 0;
	do
	{
		//oldObjValForb = objValForb;
		//oldGradForb = gradForb;

		newPar = parVect + al * dk1;
		cout << " === looking for b at\n" << newPar << "\n\n";
		calc1stOrdOptInfo( newPar, &objValForb, &gradForb );
		cout << " mult is " << gradForb.transpose() * dk1 << endl;
		al *= alStep;

		++it;
	}while( gradForb.transpose() * dk1 < 0 );

	if( gradForb.transpose() * dk1 < 0 )
	{
		cout << " something weird happened: (gradForb, dk1) < 0\n";
		return 0;
	}

	linSa = 0.0;
	//if( it > 1 && oldObjValForb <= phi0 && oldGradForb.transpose() * dk1 < 0 )
	//{
	//	cout << "~~~ better a ~~~\n";
	//	linSa = al / alStep;
	//	phiA = oldObjValForb;
	//	phiPrA = oldGradForb.transpose() * dk1;
	//}
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
		calc1stOrdOptInfo( newPar, &objValForb, &gradForb );
		phiB = objValForb;
		phiPrB = gradForb.transpose() * dk1;

		newPar = parVect + linSa * dk1;
		calc1stOrdOptInfo( newPar, &objValFora, &gradFora );
		phiA = objValFora;
		phiPrA = gradFora.transpose() * dk1;
	}
	return 1;
}


template<class PL_NUM>
void Optimizer<PL_NUM>::secant2( N_PRES* _a, N_PRES* _b )
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

template<class PL_NUM>
void Optimizer<PL_NUM>::update( N_PRES a, N_PRES b, N_PRES c, N_PRES* _a, N_PRES* _b )
{
	if( c <= a || c >= b )
	{
		*_a = a;
		*_b = b;
		return;
	}

	N_PRES objValForc;
	Matrix<N_PRES, GRAD_SIZE, 1> gradForc;
	Matrix<N_PRES, GRAD_SIZE, 1> newParc = parVect + c * dk1;

	calc1stOrdOptInfo( newParc, &objValForc, &gradForc );
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
			cout << " !!! upd while\n";
			linSd = ( 1.0 - thetaUpdate ) * aa + thetaUpdate * bb;

			N_PRES objValFord;
			Matrix<N_PRES, GRAD_SIZE, 1> gradFord;
			Matrix<N_PRES, GRAD_SIZE, 1> newPard = parVect + linSd * dk1;
			calc1stOrdOptInfo( newPard, &objValFord, &gradFord );
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

template<class PL_NUM>
void Optimizer<PL_NUM>::calc1stOrdOptInfo( const Matrix<N_PRES, GRAD_SIZE, 1>& curVal, long double* _objVal, Matrix<N_PRES, GRAD_SIZE, 1>* _gk )
{
	cout << "\tcalc 1st order\n";
	time_t begin = time( 0 );

	PL_NUM J0begin;
	PL_NUM tauBegin;
	PL_NUM B0begin;

	J0begin.elems[0] = curVal( 0 );
	J0begin.elems[1] = 1.0l;
	J0begin.elems[2] = 0.0l;
	J0begin.elems[3] = 0.0l;

	tauBegin.elems[0] = curVal( 1 );
	tauBegin.elems[1] = 0.0l;
	tauBegin.elems[2] = 1.0l;
	tauBegin.elems[3] = 0.0l;

	B0begin.elems[0] = curVal( 2 );
	B0begin.elems[1] = 0.0l;
	B0begin.elems[2] = 0.0l;
	B0begin.elems[3] = 1.0l;
		
	solver->setTask( J0begin, tauBegin, B0begin );

	cout << "\tcalculating func val\n";
	PL_NUM funcVal = calcFuncVal();
	cout << "\tfunc val done " << funcVal << endl;

	//( *_gk )( 0 ) = funcVal.elems[1] + weightJ * curVal( 0 ) / sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) );
	//( *_gk )( 1 ) = funcVal.elems[2];
	//( *_gk )( 2 ) = funcVal.elems[3] + weightJ * curVal( 2 ) / sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) );
	//*_objVal = funcVal.real() + sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) ) * weightJ;

	( *_gk )( 0 ) = funcVal.elems[1] + weightJ * curVal( 0 ) / 2.0l / ( 1.0l + M_PI * M_PI )
				* curVal( 1 ) * exp( -2.0l * charTime / curVal( 1 ) ) * ( M_PI * M_PI * ( exp( 2.0l * charTime / curVal( 1 ) ) - 1 )
				- M_PI * sin( 2.0l * M_PI * charTime / curVal( 1 ) ) + cos( 2.0l * M_PI * charTime / curVal( 1 ) ) - 1.0l );
	( *_gk )( 1 ) = funcVal.elems[2] + weightJ * curVal( 0 ) * curVal( 0 ) / 4.0l  / ( 1.0l + M_PI * M_PI ) / curVal( 1 )
				* exp( -2.0l * charTime / curVal( 1 ) ) * ( M_PI * M_PI * curVal( 1 ) * exp( 2.0l * charTime / curVal( 1 ) )
				- ( 1.0l + M_PI * M_PI ) * ( curVal( 1 ) + 2.0l * charTime ) - M_PI * curVal( 1 ) * sin( 2.0l * M_PI * charTime / curVal( 1 ) )
				+ ( curVal( 1 ) + 2.0l * ( 1.0l + M_PI * M_PI ) * charTime ) * cos( 2.0l * M_PI * charTime / curVal( 1 ) ) );
	( *_gk )( 2 ) = funcVal.elems[3] + 2.0 * weightB * curVal( 2 );
	*_objVal = funcVal.real() + weightB * curVal( 2 ) * curVal( 2 ) + weightJ * curVal( 0 ) * curVal( 0 ) / 4.0l / ( 1.0l + M_PI * M_PI )
			* curVal( 1 ) * exp( -2.0l * charTime / curVal( 1 ) ) * ( M_PI * M_PI * ( exp( 2.0l * charTime / curVal( 1 ) ) - 1 )
			- M_PI * sin( 2.0l * M_PI * charTime / curVal( 1 ) ) + cos( 2.0l * M_PI * charTime / curVal( 1 ) ) - 1.0l );

	time_t endtime = time( 0 );
	cout << "\tdone in " << endtime - begin << endl;
}

template<class PL_NUM>
void Optimizer<PL_NUM>::calc2ndOrdOptInfo( const Matrix<N_PRES, GRAD_SIZE, 1>& curVal, long double* _objVal, Matrix<N_PRES, GRAD_SIZE, 1>* _gk, Matrix<N_PRES, GRAD_SIZE, GRAD_SIZE>* _hess )
{
	cout << "\tcalc 1st order\n";
	time_t begin = time( 0 );

	PL_NUM J0begin;
	PL_NUM tauBegin;
	PL_NUM B0begin;

	J0begin.elems[0] = curVal( 0 );
	J0begin.elems[1] = 1.0l;
	J0begin.elems[2] = 0.0l;
	J0begin.elems[3] = 0.0l;

	tauBegin.elems[0] = curVal( 1 );
	tauBegin.elems[1] = 0.0l;
	tauBegin.elems[2] = 1.0l;
	tauBegin.elems[3] = 0.0l;

	B0begin.elems[0] = curVal( 2 );
	B0begin.elems[1] = 0.0l;
	B0begin.elems[2] = 0.0l;
	B0begin.elems[3] = 1.0l;
		
	solver->setTask( J0begin, tauBegin, B0begin );
	solver->calcConsts();

	cout << "\tcalculating func val\n";

	PL_NUM funcVal = calcFuncVal();

	cout << "\tfunc val done " << funcVal << endl;

	( *_gk )( 0 ) = funcVal.elems[1] + weight * curVal( 0 ) / sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) );
	( *_gk )( 1 ) = funcVal.elems[2];
	( *_gk )( 2 ) = funcVal.elems[3] + weight * curVal( 2 ) / sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) );

	*_objVal = funcVal.real() + sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) ) * weight;

	cout << " =-=- " << funcVal.real() << " " << sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) ) << " " << weight << endl;

	( *_hess )( 0, 0 ) = funcVal.elems2[0] + weight / sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) )
						* ( -curVal( 0 ) * curVal( 0 ) / ( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) ) + 1.0l );
	( *_hess )( 0, 1 ) = funcVal.elems2[1];
	( *_hess )( 0, 2 ) = funcVal.elems2[2] + weight * ( -curVal( 0 ) * curVal( 2 ) / sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) )
						/ sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) ) / sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) ) );

	( *_hess )( 1, 0 ) = funcVal.elems2[1];
	( *_hess )( 1, 1 ) = funcVal.elems2[3];
	( *_hess )( 1, 2 ) = funcVal.elems2[4];

	( *_hess )( 2, 0 ) = ( *_hess )( 0, 2 );
	( *_hess )( 2, 1 ) = funcVal.elems2[4];
	( *_hess )( 2, 2 ) = funcVal.elems2[5] + weight / sqrt( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) )
						* ( -curVal( 2 ) * curVal( 2 ) / ( curVal( 0 ) * curVal( 0 ) + curVal( 2 ) * curVal( 2 ) ) + 1.0l );

	time_t endtime = time( 0 );
	cout << "\tdone in " << endtime - begin << endl;
}

template<class PL_NUM>
PL_NUM Optimizer<PL_NUM>::calcFuncVal()
{
	//PL_NUM funcVal( 0.0, 0.0 );
	PL_NUM funcVal;

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

template<class PL_NUM>
N_PRES Optimizer<PL_NUM>::calcBettaN()
{
	Matrix<N_PRES, GRAD_SIZE, 1> yk = gk1 - gk;

	return 1.0 / ( dk.transpose() * yk ) * ( yk - 2.0 * dk * ( yk.transpose() * yk ) / ( dk.transpose() * yk ) ).transpose() * gk1;
}

template<class PL_NUM>
N_PRES Optimizer<PL_NUM>::calcBettaN_()
{
	N_PRES oldBetta = calcBettaN();
	N_PRES etta = 0.01;		//as in Hager's paper
	N_PRES ettak = -1.0 / ( dk.norm() * min( etta, gk.norm() ) );
	
	return max( oldBetta, ettak );
}

template<class PL_NUM>
N_PRES Optimizer<PL_NUM>::min( N_PRES a, N_PRES b )
{
	N_PRES ret = a;
	if( b < ret )
	{
		ret = b;
	}
	return ret;
}

template<class PL_NUM>
N_PRES Optimizer<PL_NUM>::max( N_PRES a, N_PRES b )
{
	N_PRES ret = a;
	if( b > ret )
	{
		ret = b;
	}
	return ret;
}

#endif