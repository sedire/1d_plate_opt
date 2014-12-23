#include "ParamSack.h"

N_PRES ParamSack::getE1()
{
	return E1;
}
N_PRES ParamSack::getE2()
{
	return E2;
}
N_PRES ParamSack::getNu21()
{
	return nu21;
}
N_PRES ParamSack::getRho()
{
	return rho;
}
N_PRES ParamSack::getSigma_x()
{
	return sigma_x;
}
N_PRES ParamSack::getMu()
{ 
	return mu;
}
N_PRES ParamSack::getEps_0()
{
	return eps_0;
}
N_PRES ParamSack::getEps_x()
{
	return eps_x;
}
N_PRES ParamSack::getH()
{
	return h;
}
N_PRES ParamSack::getA()
{
	return a;
}
int ParamSack::getNumberOfScenarios()
{
	return numberOfScenarios;
}
int ParamSack::getCurrentType()
{
	return currentType;
}
int ParamSack::getNumberOfCurrentParams()
{
	return numberOfCurrentParams;
}
N_PRES ParamSack::getCurrentParams( int scen, int paramNum )
{
	return currentParams[scen][paramNum];
}
int ParamSack::getStressType()
{
	return stressType;
}
int ParamSack::getNumberOfStressParams()
{
	return numberOfStressParams;
}
N_PRES ParamSack::getStressParams( int scen, int paramNum )
{
	return stressParams[scen][paramNum];
}
N_PRES ParamSack::getTotalTimeT1()
{
	return totalTimeT1;
}
N_PRES ParamSack::getSwitchTimeT0()
{
	return switchTimeT0;
}
int ParamSack::getNodeNumOnY()
{
	return nodeNumOnY;
}
N_PRES ParamSack::getDt()
{
	return dt;
}
int ParamSack::getMaxNewtonIter()
{
	return maxNewtonIter;
}
N_PRES ParamSack::getOrthonormCheckEps()
{
	return orthonormCheckEps;
}
N_PRES ParamSack::getNewtonEps()
{
	return newtonEps;
}
N_PRES ParamSack::getNewtonAlmostZero()
{
	return newtonAlmostZero;
}

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
		ss >> word;
		if( word.compare( "E1" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> E1;
			cout << E1 << endl;
		}

		getline( ifs, line );
		ss.clear();	//clear any bits set
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "E2" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> E2;
			cout << E2 << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "nu21" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> nu21;
			cout << nu21 << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "rho" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> rho;
			cout << rho << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "sigma_x" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> sigma_x;
			cout << sigma_x << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "mu" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> mu;
			cout << mu << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "eps_0" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> eps_0;
			cout << eps_0 << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "eps_x" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> eps_x;
			cout << eps_x << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "h" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> h;
			cout << h << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "a" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> a;
			cout << a << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "numberOfScenarios" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> numberOfScenarios;
			cout << numberOfScenarios << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "currentType" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			string currentTypeInput;
			ss >> currentTypeInput;
			cout << currentTypeInput << endl;
			if( currentTypeInput.compare( "currentConst" ) == 0 )
			{
				currentType = currentConst;
			}
			else if( currentTypeInput.compare( "currentSin" ) == 0 )
			{
				currentType = currentSin;
			}
			else if( currentTypeInput.compare( "currentExpSin" ) == 0 )
			{
				currentType = currentExpSin;
			}
			else
			{
				cout << "ERROR: file " << inputFname << " has bad format\n";
				return;
			}
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "numberOfCurrentParams" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> numberOfCurrentParams;
			cout << numberOfCurrentParams << endl;
			currentParams.resize( numberOfScenarios + 1, vector<N_PRES>( numberOfCurrentParams, 0.0 ) );
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "currentParamsFirstStage" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			for( int i = 0; i < numberOfCurrentParams; ++i )
			{
				ss >> currentParams[0][i];
				cout << currentParams[0][i] << endl;
			}
		}

		for( int scen = 0; scen < numberOfScenarios; ++scen )
		{
			getline( ifs, line );
			ss.clear();
			ss.str( std::string() );
			ss << line;
			ss >> word;
			if( word.compare( "currentParamsSecondStage" ) != 0 )
			{
				cout << "ERROR: file " << inputFname << " has bad format\n";
				return;
			}
			else
			{
				for( int i = 0; i < numberOfCurrentParams; ++i )
				{
					ss >> currentParams[scen][i];
					cout << currentParams[scen][i] << endl;
				}
			}
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "stressType" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			string stressTypeInput;
			ss >> stressTypeInput;
			cout << stressTypeInput << endl;
			if( stressTypeInput.compare( "stressWhole" ) == 0 )
			{
				stressType = stressWhole;
			}
			else if( stressTypeInput.compare( "stressCentered" ) == 0 )
			{
				stressType = stressCentered;
			}
			else
			{
				cout << "ERROR: file " << inputFname << " has bad format\n";
				return;
			}
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "numberOfStressParams" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> numberOfStressParams;
			cout << numberOfStressParams << endl;
			stressParams.resize( numberOfScenarios, vector<N_PRES>( numberOfStressParams, 0.0 ) );
		}

		for( int scen = 0; scen < numberOfScenarios; ++scen )
		{
			getline( ifs, line );
			ss.clear();
			ss.str( std::string() );
			ss << line;
			ss >> word;
			if( word.compare( "stressParams" ) != 0 )
			{
				cout << "ERROR: file " << inputFname << " has bad format\n";
				return;
			}
			else
			{
				for( int i = 0; i < numberOfStressParams; ++i )
				{
					ss >> stressParams[scen][i];
					cout << stressParams[scen][i] << endl;
				}
			}
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "totalTimeT1" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> totalTimeT1;
			cout << totalTimeT1 << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "switchTimeT0" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> switchTimeT0;
			cout << switchTimeT0 << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "nodeNumOnY" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> nodeNumOnY;
			cout << nodeNumOnY << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "dt" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> dt;
			cout << dt << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "maxNewtonIter" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> maxNewtonIter;
			cout << maxNewtonIter << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "orthonormCheckEps" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> orthonormCheckEps;
			cout << orthonormCheckEps << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "newtonEps" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> newtonEps;
			cout << newtonEps << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "newtonAlmostZero" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return;
		}
		else
		{
			ss >> newtonAlmostZero;
			cout << newtonAlmostZero << endl;
		}
	}
	else
	{
		cout << "ERROR: cannot open file " << inputFname << " to load parameter into ParamSack\n";
		return;
	}
	return;
}