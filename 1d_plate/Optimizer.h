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
#include "ParamSack.h"

using namespace Eigen;
using std::cout;
using std::endl;
using std::ofstream;

extern ParamSack GlobalParams;

//double calc1stOrdOptInfoCG_DES( double* g, double* x, long n );

void optimizeASA_Taus();
void optimizeASAPiece();

#endif