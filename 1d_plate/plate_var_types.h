#ifndef _PLATE_1D_PLATE_VAR_TYPES_
#define _PLATE_1D_PLATE_VAR_TYPES_ 1

#define _USE_MATH_DEFINES

#define THREAD_NUM 2

#define N_PRES long double
#define EQ_NUM 8
//#define NEWTON_ALMOST_ZERO 1e-11l
//#define NEWTON_EPS 1e-6l
#define J0_SCALE 100000000.0l
#define CHAR_TIME 0.05l
#define SWITCH_TIME 0.01l
#define DELTA_T 0.0001
#define NODES_Y 10001	//CHANGE THIS BACK to 10001
#define ABM_STAGE_NUM 4
#define W_SCALE 1.0l

#define ORTHONORM_CHECK_EPS 1e-7

#define J_WEIGHT 1.0l

#define GRAD_SIZE 3
#define GRAD_SIZE_FULL 12
#define GRAD_SIZE_FIRST 3
#define GRAD_SIZE_SECOND 6

#define MAX_NEWTON_ITER 10

#define SCEN_NUMBER 3

enum {stressWhole, stressCentered};
enum {currentConst, currentSin, currentExpSin, currentPieceLin};

const N_PRES GlobalP01( 750000 );
const N_PRES GlobalP02( 1000000 );
const N_PRES GlobalP03( 2000000 );

const N_PRES GlobalTauP1( 0.008 );
const N_PRES GlobalTauP2( 0.01 );
const N_PRES GlobalTauP3( 0.012 );

//global arrays for adjoint
extern N_PRES* GlobalResArrays;
extern N_PRES* GlobalResDtArrays;
extern N_PRES* GlobalResAdjArrays;

#endif