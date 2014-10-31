#ifndef _PLATE_1D_HAGEROPTFUNCS_
#define _PLATE_1D_HAGEROPTFUNCS_ 1

#include "hyperDual.h"
#include "time.h"
#include "Plate.h"
#include "Solver.h"
#include "plate_var_types.h"
#include <iostream>

#include "asa_user.h"
#include "AdjSolver.h"

using std::cout;
using std::endl;

double calc1stOrdOptInfoCG_DES( double* g, double* x, long n );
double calcValGradTaus( double* g, double* x, long n );
double calcValGradTausAdj( double* g, double* x, long n );			//for the objective function splitted into 2 integrals with T-based weights
double calcValGradTausAdjSolid( double* g, double* x, long n );		//for the objective function that is not splitted into 2 integrals with T-based weights
double calcValTaus( double* x, long n );

//double calcValCG_DES( double* x, long n );

//void calcGradCG_DES( double* g, double* x, long n );

double calc1stOrdOptInfoASA( asa_objective* asa );
double calcValASA( asa_objective* asa );
void calcGradASA( asa_objective* asa );
///////////////////////////////////////
double calc1stOrdOptInfoASA_Taus( asa_objective* asa );
double calcValASA_Taus( asa_objective* asa );
void calcGradASA_Taus( asa_objective* asa );



#endif