#ifndef __CONGRADMIN_H__
#define __CONGRADMIN_H__
   #include "tzmatlab.h"

   #include <string.h>




   void frprmn(double p[], int n, int *iter, double *fret,
               double (*func)(double [], int), void (*dfunc)(double [], double [], int, double (*func)(double [], int), double *, double),
               double *ftol_p, int *itmax_p, double *tol_brent_p, int *itmax_brent_p, double *grdh_p);
/*        //Outputs:   ansi-c*/
/*        //  p[0, ..., n-1]:  the location of the minimum if it converges, which replaces the starting value.   ansi-c*/
/*        //  iter:  pointer to the number of iterations that were performed.   ansi-c*/
/*        //  fret:  pointer to the minimum value of the function.   ansi-c*/
/*        //Inputs:   ansi-c*/
/*        //  p[0, ..., n-1]:  a starting point for the minimization.   ansi-c*/
/*        //  n:  the dimension of p.   ansi-c*/
/*        //  ftol_p:  pointer to the convergence tolerance on the objective function value. Default: 1.0e-4 if NULL.   ansi-c*/
/*        //  itmax_p:    pointer to the maximum number of iterations in the main minimization program frprmn().  Default: 2000 if NULL.   ansi-c*/
/*        //  tol_brent_p:  pointer to the convergence tolerance for the line minimization in brent().  Default: 2.0e-4 if NULL.   ansi-c*/
/*        //  itmax_brent_p:  pointer to the maximum number of iterations for the line minimization in brent().  Default: 100 if NULL.   ansi-c*/
/*        //  grdh:  pointer to the user's specified step size for a numerical gradient.  If NULL, dfunc() (i.e., gradcd_gen()) will select grdh automatically.   ansi-c*/
/*        //  func():  the objective function.   ansi-c*/
/*        //  dfunc(): the gradient function computing the numerical gradient.  In the form of gradcd_gen() in cstz.c.   ansi-c*/

   void congradmin_SetPrintFile(char *filename);
/*        //If filename=NULL, no intermediate results will be printed out to a file.   ansi-c*/
/*  //   void congradmin_SetPrintFile(FILE *fptr_sp);   ansi-c*/
/*        //If fptr_sp=NULL, no intermediate results will be printed out to a file.   ansi-c*/
/*  //   void congradmin_SetPrintFile_db(FILE *fptr_sp);   ansi-c*/
#endif
