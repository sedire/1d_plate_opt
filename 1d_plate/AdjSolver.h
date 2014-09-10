#include "plate_var_types.h" 
#include <iostream>

using std::cout;
using std::endl;


class AdjSolver		
{
public:
	AdjSolver();
	~AdjSolver();
	N_PRES do_step();		//the adjoint problem is solved backwards in time
	void decreaseTime();

private:
	N_PRES dt;
	N_PRES curTime;
	int curTimeStep;
};

void AdjSolver::decreaseTime()
{
	--curTimeStep;
	curTime -= dt;
}