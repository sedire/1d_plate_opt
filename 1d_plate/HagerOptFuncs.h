#ifndef _PLATE_1D_HAGEROPTFUNCS_
#define _PLATE_1D_HAGEROPTFUNCS_ 1

#include "asa_user.h"

double calc1stOrdOptInfoCG_DES( double* g, double* x, long n );
double calcValGradTaus( double* g, double* x, long n );
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