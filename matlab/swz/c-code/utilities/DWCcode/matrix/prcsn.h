
/*
   Defines the precision to be used
*/

#ifndef __PRECISION_H__
#define __PRECISION_H__

#include <float.h>

/********** double precision **********/
#define PRECISION              double
#define MACHINE_EPSILON        1.11E-16
#define SQRT_MACHINE_EPSILON   1.06E-08
#define PRECISION_SIZE         8
#define PRECISION_SHIFT        3
#define PRECISION_WORD         qword
#define MINUS_INFINITY        -1.0E300
#define PLUS_INFINITY          1.0E300
//#define MINUS_INFINITY       -DBL_MAX
//#define PLUS_INFINITY         DBL_MAX
/**************************************/

/********** single precision **********
#define PRECISION             float
#define MACHINE_EPSILON       5.97E-08
#define SQRT_MACHINE_EPSILON  2.45E-04
#define PRECISION_SIZE        4
#define PRECISION_SHIFT       2
#define PRECISION_WORD        dword
#define MINUS_INFINITY       -FLT_MAX
#define PLUS_INFINITY         FLT_MAX
/**************************************/

#endif
