#ifndef _PLATE_1D_PLATE_VAR_TYPES_
#define _PLATE_1D_PLATE_VAR_TYPES_ 1

#define _USE_MATH_DEFINES

#define THREAD_NUM 4

#define N_PRES long double
//#define PL_NUM complex<long double>
#define EQ_NUM 8
#define ALMOST_ZERO 0.0000000001l
#define J0_SCALE 100000000.0l
#define BY0_SCALE 1.0l
#define CHAR_TIME 0.05l
#define SWITCH_TIME 0.01l
#define DELTA_T 0.0001
#define NODES_Y 10001

#define J_WEIGHT 1.0l

#define GRAD_SIZE 3
#define GRAD_SIZE_FULL 12
#define GRAD_SIZE_FIRST 3
#define GRAD_SIZE_SECOND 6

#define MAX_NEWTON_ITER 10

#define SCEN_NUMBER 3

enum {stress_whole, stress_centered};
enum {current_const, current_sin, current_exp_sin};

const N_PRES GlobalP01( 7500000 );
const N_PRES GlobalP02( 10000000 );
const N_PRES GlobalP03( 20000000 );

const N_PRES GlobalTauP1( 0.008 );
const N_PRES GlobalTauP2( 0.01 );
const N_PRES GlobalTauP3( 0.012 );

#endif