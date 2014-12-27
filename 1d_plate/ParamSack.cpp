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
N_PRES ParamSack::getBy0()
{
	return By0;
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
N_PRES ParamSack::getCurrentParams1st( int paramNum )
{
	return currentParams1st[paramNum];
}
N_PRES ParamSack::getCurrentParams2nd( int scen, int paramNum )
{
	return currentParams2nd[scen][paramNum];
}
vector<N_PRES> ParamSack::getCurrentParams( int scen )
{
	vector<N_PRES> ret;
	ret.resize( numberOfCurrentParams, 0.0 );
	for( int i = 0; i < currentParams1st.size(); ++i )
	{
		ret[i] = currentParams1st[i];
	}
	for( int i = 0; i < currentParams2nd[scen].size(); ++i )
	{
		ret[currentParams1st.size() + i] = currentParams2nd[scen][i];
	}
	return ret;
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
vector<N_PRES> ParamSack::getStressParams( int scen )
{
	vector<N_PRES> ret;
	ret.resize( numberOfStressParams, 0.0 );
	for( int i = 0; i < numberOfStressParams; ++i )
	{
		ret[i] = stressParams[scen][i];
	}
	return ret;
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

	By0( -1.0 ),

	currentType( -1 ),
	numberOfCurrentParams( -1 ),

	stressType( -1 ),
	numberOfStressParams( -1 ),

	totalTimeT1( 0.0 ),
	switchTimeT0( 0.0 ),

	//parameters if the numerical methods
	nodeNumOnY( -1 ),
	dt( 0.0 ),
	maxNewtonIter( -1 ),
	orthonormCheckEps( 0.0 ),
	newtonEps( 0.0 ),
	newtonAlmostZero( 0.0 )
{
}

ParamSack::~ParamSack()
{
}

int ParamSack::loadFromFile( string inputFname )
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
			return 1;
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
			return 1;
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
			return 1;
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
			return 1;
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
			return 1;
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
			return 1;
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
			return 1;
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
			return 1;
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
			return 1;
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
			return 1;
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
		if( word.compare( "totalTimeT1" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return 1;
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
			return 1;
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
		if( word.compare( "By0" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return 1;
		}
		else
		{
			ss >> By0;
			cout << By0 << endl;
		}

		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "numberOfScenarios" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return 1;
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
			return 1;
		}
		else
		{
			string currentTypeInput;
			ss >> currentTypeInput;
			cout << currentTypeInput << endl;
			if( currentTypeInput.compare( "currentExpSin" ) == 0 )
			{
				currentType = currentExpSin;
			}
			else if( currentTypeInput.compare( "currentPieceLin" ) == 0 )
			{
				currentType = currentPieceLin;
			}
			else
			{
				cout << "ERROR: file " << inputFname << " has bad format\n";
				return 1;
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
			return 1;
		}
		else
		{
			ss >> numberOfCurrentParams;
			cout << numberOfCurrentParams << endl;
			if( currentType == currentExpSin && numberOfCurrentParams != 6 )
			{
				cout << "ERROR: wrong number of parameters for currentExpSin\n";
				return 1;
			}
			else if( currentType == currentPieceLin && numberOfCurrentParams < 6 )
			{
				cout << "ERROR: wrong number of parameters for currentPieceLine\n";
				return 1;
			}
		}

		int numberOfCurrentParams1st = -1;
		getline( ifs, line );
		ss.clear();
		ss.str( std::string() );
		ss << line;
		ss >> word;
		if( word.compare( "currentParamsFirstStage" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return 1;
		}
		else
		{
			if( currentType == currentExpSin )
			{
				numberOfCurrentParams1st = 3;
			}
			else if( currentType == currentPieceLin )
			{
				N_PRES dT = totalTimeT1 / numberOfCurrentParams;
				numberOfCurrentParams1st = (int)floor( switchTimeT0 / dT );
			}
			else
			{
				cout << "ERROR: this type of current is not supported yet\n";
				return 1;
			}
			currentParams1st.resize( numberOfCurrentParams1st, 0.0 );
			for( int i = 0; i < numberOfCurrentParams1st; ++i )
			{
				ss >> currentParams1st[i];
				cout << currentParams1st[i] << endl;
			}
		}

		for( int scen = 0; scen < numberOfScenarios; ++scen )
		{
			int numberOfCurrentParams2nd = numberOfCurrentParams - numberOfCurrentParams1st;
			currentParams2nd.resize( numberOfScenarios, vector<N_PRES>( numberOfCurrentParams2nd, 0.0 ) );
			getline( ifs, line );
			ss.clear();
			ss.str( std::string() );
			ss << line;
			ss >> word;
			if( word.compare( "currentParamsSecondStage" ) != 0 )
			{
				cout << "ERROR: file " << inputFname << " has bad format\n";
				return 1;
			}
			else
			{
				for( int i = 0; i < numberOfCurrentParams2nd; ++i )
				{
					ss >> currentParams2nd[scen][i];
					cout << currentParams2nd[scen][i] << endl;
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
			return 1;
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
				return 1;
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
			return 1;
		}
		else
		{
			ss >> numberOfStressParams;
			cout << numberOfStressParams << endl;
			if( stressType == stressWhole && numberOfStressParams != 1 )
			{
				cout << "ERROR: file " << inputFname << " wrong number of stress params\n";
				return 1;
			}
			else if( stressType == stressCentered && numberOfStressParams != 3 )
			{
				cout << "ERROR: file " << inputFname << " wrong number of stress params\n";
				return 1;
			}
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
				return 1;
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
		if( word.compare( "nodeNumOnY" ) != 0 )
		{
			cout << "ERROR: file " << inputFname << " has bad format\n";
			return 1;
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
			return 1;
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
			return 1;
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
			return 1;
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
			return 1;
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
			return 1;
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
		return 1;
	}
	return 0;
}