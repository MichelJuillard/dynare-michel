#ifndef __CONGRADMIN_H__
#define __CONGRADMIN_H__
   #include "tzmatlab.h"

   #include <string.h>




   void frprmn(double p[], int n, int *iter, double *fret,
               double (*func)(double [], int), void (*dfunc)(double [], double [], int, double (*func)(double [], int), double *, double),
               double *ftol_p, int *itmax_p, double *tol_brent_p, int *itmax_brent_p, double *grdh_p);
      //Outputs:
      //  p[0, ..., n-1]:  the location of the minimum if it converges, which replaces the starting value.
      //  iter:  pointer to the number of iterations that were performed.
      //  fret:  pointer to the minimum value of the function.
      //Inputs:
      //  p[0, ..., n-1]:  a starting point for the minimization.
      //  n:  the dimension of p.
      //  ftol_p:  pointer to the convergence tolerance on the objective function value. Default: 1.0e-4 if NULL.
      //  itmax_p:    pointer to the maximum number of iterations in the main minimization program frprmn().  Default: 2000 if NULL.
      //  tol_brent_p:  pointer to the convergence tolerance for the line minimization in brent().  Default: 2.0e-4 if NULL.
      //  itmax_brent_p:  pointer to the maximum number of iterations for the line minimization in brent().  Default: 100 if NULL.
      //  grdh:  pointer to the user's specified step size for a numerical gradient.  If NULL, dfunc() (i.e., gradcd_gen()) will select grdh automatically.
      //  func():  the objective function.
      //  dfunc(): the gradient function computing the numerical gradient.  In the form of gradcd_gen() in cstz.c.

   void congradmin_SetPrintFile(char *filename);
      //If filename=NULL, no intermediate results will be printed out to a file.
//   void congradmin_SetPrintFile(FILE *fptr_sp);
      //If fptr_sp=NULL, no intermediate results will be printed out to a file.
//   void congradmin_SetPrintFile_db(FILE *fptr_sp);
#endif
