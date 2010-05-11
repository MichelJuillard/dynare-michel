#ifndef __CSMINWEL_H__
#define __CSMINWEL_H__

#include "tzmatlab.h"

#include <string.h>
#include <float.h>

/*  //--- This extern variable allows an input by the user from an input data file.   ansi-c*/
extern double GRADSTPS_CSMINWEL;

void csminwel(double (*fcn)(double *x, int n, double **args, int *dims),
            double *x, int n, double *H, double *gh,
            int (*grad)(double *x, int n, double *g, double **args, int *dims),
            double *fh, double crit, int *itct, int nit,
            int *fcount, int *retcodeh, double **args, int *dims);
/*  // Alternative but less clear way:  ... (double (*fcn)(double *, int, double **, int *), ...   ansi-c*/

void csminwel_SetPrintFile(char *filename);
int csminwel_randomseedChanged(int seednumber);


#endif
