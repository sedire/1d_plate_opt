#include "ParamSack.h"

ParamSack::ParamSack() :
	E1( 0.0 ),
	E2( 0.0 ),
	nu21( 0.0 ),			//Poisson's ratio	
	rho( 0.0 ),				//composite's density
	sigma_x( 0.0 ),			//electric conductivity
	mu( 0.0 ),
	eps_0( 0.0 ),
	eps_x( 0.0 ),

	h( 0.0 ),				//thickness of the plate
	a( 0.0 ),				//width of the plate

	currentType( -1 ),
	numberOfCurrentParams( 0 ),

	stressType( -1 ),
	numberOfStressParams( 0 ),

	totalTimeT1( 0.0 ),
	switchTimeT0( 0.0 ),

	//parameters if the numerical methods
	nodeNumOnY( 0 ),
	dt( 0.0 ),
	maxNewtonIter( 0 ),
	orthonormCheckEps( 0.0 ),
	newtonEps( 0.0 ),
	newtonAlmostZero( 0.0 )
{
}

ParamSack::~ParamSack()
{
}

void ParamSack::loadFromFile( string inputFname )
{
	ifstream ifs( inputFname );
	if( ifs.is_open() )
	{
		string line;

		getline( ifs, line );
		stringstream ss( line );
		string word;
		getline( ss, word, ' ' );
		if( word.compare( "E1" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			getline( ss, word, ' ' );
			cout << word << endl;
		}
	}
	else
	{
		cout << "ERROR: cannot open file " << inputFname << " to load parameter into ParamSack\n";
		return;
	}
	return;
}