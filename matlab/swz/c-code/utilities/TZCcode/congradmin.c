/*************************************************************
 *  Conjugate Gradient Minimization Methods.  See Numerical Recipes in C by Press, Flannery, Teukolsky, and Vetterling.
 *  (I)  frprmn():  Plolak-Ribiere method with the line minimization without using the derivative information.
 *  (II) dlinmin():  Fletcher-Reeves method with the line minimization using the derivative information.
 *
 * Modified by Tao Zha, September 2003.
*************************************************************/

#include "congradmin.h"

static void linmin(double p[], double xi[], int n, double *fret, double tol_brent, int itmax_brent, double (*func)(double [], int));
static double brent(double ax, double bx, double cx, double (*f)(double), double tol_brent, double itmax_brent, double *xmin);
//
static void dlinmin(double p[], double xi[], int n, double *fret, double tol_dbrent, double itmax_dbrent, double *grdh_p, double (*func)(double [], int), void (*dfunc)(double [], double [], int, double (*func)(double [], int), double *, double));
static double dbrent(double ax, double bx, double cx, double (*f)(double), double (*df)(double, double *), double *grdh_p, double tol_dbrent, double itmax_dbrent, double *xmin);
static double df1dim(double x, double *grdh_p);
//
static void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double));
static double f1dim(double x);
//
static double ftd_norm2(double *vnew_p, double *vold_p, int _n);
static double ftd_innerproduct(double *x, double *y, int _n);


#define ANGLE 0.001  //.0   implies 90.00 degress (acrcos(ANGLE)*180/pi).
                     //.005 implies 89.71 degrees (acrcos(ANGLE)*180/pi).
                     //.01  implies 89.43 degrees (acrcos(ANGLE)*180/pi).
                     //.05  implies 87.13 degrees (acrcos(ANGLE)*180/pi).
                     //.1   implies 84.26 degrees (acrcos(ANGLE)*180/pi).
#define STRLEN 192
static FILE *fptr_interesults = (FILE *)NULL;   //Printing intermediate results to a file.
static char filename_sp3vecs[STRLEN];  //Three vectors.  1st row: line search direction; 2nd row: numerical gradient; 3rd row: vectorized parameters.
//static FILE *fptr_interesults_db = (FILE *)NULL;   //Printing intermediate results to a file for debugging (db).
#define PRINTON    //Added by TZ, September 2003.
#define EPS 1.0e-10        //Small number to rectify special case of converging to exactly zero function value.
#ifdef PRINTON      //Added by TZ, September 2003.
   #define FREEALL {tzDestroy(xi); tzDestroy(h); tzDestroy(g); tzDestroy(pold); tzDestroy(numgrad)}
#else
   #define FREEALL {tzDestroy(xi); tzDestroy(h); tzDestroy(g);}
#endif
void frprmn(double p[], int n, int *iter, double *fret,
            double (*func)(double [], int), void (*dfunc)(double [], double [], int, double (*func)(double [], int), double *, double),
			double *ftol_p, int *itmax_p, double *tol_brent_p, int *itmax_brent_p, double *grdh_p) {
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
   int j, its, itmax, itmax_brent;
   double gg, gam, fp, dgg, ftol, tol_brent;
   double *g=NULL, *h=NULL, *xi=NULL;
   #ifdef PRINTON      //Added by TZ, September 2003.
   time_t begtime, currentime;
   double normforp, *pold = NULL, *numgrad = NULL;
   int cnt_wrong_dirs = -1;  //Counts the number of times that a numerical direction in the line search has a wrong sign.
   #endif

   //=== Memory allocation.
   g=tzMalloc(n, double);
   h=tzMalloc(n, double);
   xi=tzMalloc(n, double);
   //
   numgrad = tzMalloc(n, double);   //Added by TZ, September 2003.
   #ifdef PRINTON      //Added by TZ, September 2003.
      pold = tzMalloc(n, double);
   #endif

   //=== Default values.
   if (!ftol_p)  ftol = 1.0e-4;  else  ftol = *ftol_p;
   if (!itmax_p)  itmax = 200;  else  itmax = *itmax_p;
   if (!tol_brent_p)  tol_brent = 2.0e-4;  else  tol_brent = *tol_brent_p;
   if (!itmax_brent_p)  itmax_brent = 100;  else  itmax_brent = *itmax_brent_p;

   fp=(*func)(p, n);
   (*dfunc)(xi, p, n, func, grdh_p, fp);
   for (j=n-1;j>=0;j--) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
   memcpy(numgrad, xi, n*sizeof(double));   //Added by TZ, September 2003. Save the numerical gradient to be printed out at the right place.
   for (its=0;its<itmax;its++) {
      #ifdef PRINTON
      time(&begtime);    //Added by TZ, September 2003.
      memcpy(pold, p, n*sizeof(double));   //Added by TZ, September 2003.
      #endif
      //====== Added by TZ, September 2003 ======
      if ( !(fptr_interesults = fopen(filename_sp3vecs,"w")) ) {
         printf("\n\nUnable to create the starting point data file %s in congradmin.c!\n", filename_sp3vecs);
         getchar();
         exit(EXIT_FAILURE);
      }
      // rewind(fptr_interesults);   //Must put the pointer at the beginning of the file.
      //=== Prints out the line search direction.
      fprintf(fptr_interesults, "--------Line search direction---------\n");
      for (j=0; j<n; j++)  fprintf(fptr_interesults, " %0.16e ", xi[j]);
      fprintf(fptr_interesults, "\n");
      // fflush( fptr_interesults );
      //=== Prints out the message about a wrong numerical direction in the line search for the miminziation.
      if ( ftd_innerproduct(xi, numgrad, n)/(ftd_norm2(xi, xi, n)*ftd_norm2(numgrad, numgrad, n)) > - ANGLE ) {
         #ifdef PRINTON
         printf("\n----------------\n"
                 "Warning: wrong numerical direction in the line search for the miminziation (a total of %d times)!\n"
                 "----------------\n", ++cnt_wrong_dirs);
         #endif
      }


      *iter=its;
      #if defined (CGI_OPTIMIZATION)
         linmin(p,xi,n,fret, tol_brent, itmax_brent, func);
      #elif defined (CGII_OPTIMIZATION)
         dlinmin(p, xi, n, fret, tol_brent, itmax_brent, grdh_p, func, dfunc);
      #else
         fn_DisplayError("The minimization routine frprmn() requires activating CGI_OPTIMIZATION or CGII_OPTIMIZATION in tzmatlab.h");
      #endif
      #ifdef PRINTON
         normforp = ftd_norm2(p, pold, n);
         //=== Prints out intermediate results.
         printf("\n========================================\n");
         printf("Intermediate results for the conjugate gradient algorithm.");
         printf("\n (1) Number of iterations so far (maximum number): %d (%d)\n (2) New value of objective function (old value, improvement): %0.9f (%0.9f, %0.9f)\n"
                " (3) Norm-2 of dx: %0.9f\n",
                its, itmax, *fret, fp, fp-(*fret), normforp);
         fflush(stdout);                // Flush the buffer to get out this message without delay.
      #endif
      //====== The following statements print out intermediate results.  Added by TZ, September 2003 ======
      //=== Prints out the gradient.
      fprintf(fptr_interesults, "--------Numerical graident---------\n");
      for (j=0; j<n; j++)  fprintf(fptr_interesults, " %0.16e ", numgrad[j]);
      fprintf(fptr_interesults, "\n");
      //
      fprintf(fptr_interesults, "--------Restarting point---------\n");
      for (j=0; j<n; j++)  fprintf(fptr_interesults, " %0.16e ", p[j]);
      fprintf(fptr_interesults, "\n\n");
//      fflush( fptr_interesults );
      tzFclose(fptr_interesults);


      if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
         //This is a normal convergence.
         printf("\n----- Normal convergence by the criterion of the objective function evaluation -----------\n");
			FREEALL
			return;
		}
      fp=(*func)(p, n);
      (*dfunc)(xi, p, n, func, grdh_p, fp);
      memcpy(numgrad, xi, n*sizeof(double));   //Added by TZ, September 2003. Save the numerical gradient to be printed out at the right place.
//      if (filename_sp3vecs) {
//         //=== Prints out the gradient.
//         fprintf(fptr_interesults, "--------Numerical graident---------\n");
//         for (j=0; j<n; j++)  fprintf(fptr_interesults, " %0.16e ", xi[j]);
//         fprintf(fptr_interesults, "\n\n");
////         fflush( fptr_interesults );

//         tzFclose(fptr_interesults);
//      }
      dgg=gg=0.0;
      for (j=n-1;j>=0;j--) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
      for (j=n-1;j>=0;j--) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}

      #ifdef PRINTON
      time(&currentime);
      //=== Times the iterative progress.
      printf(" (4) Seconds to complete one iteration: %0.4f\n (5) Current time of day: %s\n", difftime(currentime, begtime), ctime(&currentime));
      fflush(stdout);                // Flush the buffer to get out this message without delay.
      #endif
	}
   fn_DisplayError("The maximum number of iterations in frprmn() is reached before convergence");
}
#undef PRINTON
#undef EPS
#undef FREEALL


#if defined (CGI_OPTIMIZATION)
   static int ncom;
   static double *pcom=NULL, *xicom=NULL, (*nrfunc)(double [], int);   //nrfunc(), pcom, ncom, and xicom will be used by f1dim().
   static void linmin(double p[], double xi[], int n, double *fret, double tol_brent, int itmax_brent, double (*func)(double [], int)) {
      //Outputs:
      //  p[0, ..., n-1]:  a returned and reset value.
      //  xi[0, ..., n-1]:  a value repaced by the actual vector displacement that p was moved.
      //  fret:  the value of func at the returned location p.
      //Inputs:
      //  p[0, ..., n-1]:  a given point.
      //  xi[0, ..., n-1]:  a given multidimensional direction.
      //  n:  the dimension of p and xi.
      //  func():  the objective function.
      int j;
      double xx,xmin,fx,fb,fa,bx,ax;

      ncom=n;
      pcom = tzMalloc(n, double);
      xicom = tzMalloc(n, double);
      nrfunc=func;
      for (j=n-1;j>=0;j--) {
         pcom[j]=p[j];
         xicom[j]=xi[j];
      }
      ax=0.0;
      xx=1.0;
      mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
      *fret=brent(ax,xx,bx,f1dim, tol_brent, itmax_brent, &xmin);
      for (j=n-1;j>=0;j--) {
         xi[j] *= xmin;
         p[j] += xi[j];
      }
      tzDestroy(xicom);
      tzDestroy(pcom);
   }


   //=== Used by linmin() only;
   #define CGOLD 0.3819660
   #define ZEPS 1.0e-10
   #define SHFT(a,b,c,d)  {(a)=(b);(b)=(c);(c)=(d);}
   #define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
   static double brent(double ax, double bx, double cx, double (*f)(double), double tol_brent, double itmax_brent, double *xmin) {
      int iter;
      double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
      double e=0.0;

      a=(ax < cx ? ax : cx);
      b=(ax > cx ? ax : cx);
      x=w=v=bx;
      fw=fv=fx=(*f)(x);
      for (iter=0;iter<itmax_brent;iter++) {
         xm=0.5*(a+b);
         tol2=2.0*(tol1=tol_brent*fabs(x)+ZEPS);
         if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
            *xmin=x;
            return fx;
         }
         if (fabs(e) > tol1) {
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if (q > 0.0) p = -p;
            q=fabs(q);
            etemp=e;
            e=d;
            if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
               d=CGOLD*(e=(x >= xm ? a-x : b-x));
            else {
               d=p/q;
               u=x+d;
               if (u-a < tol2 || b-u < tol2)
                  d=SIGN(tol1,xm-x);
            }
         } else {
            d=CGOLD*(e=(x >= xm ? a-x : b-x));
         }
         u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
         fu=(*f)(u);
         if (fu <= fx) {
            if (u >= x) a=x; else b=x;
            SHFT(v,w,x,u)
            SHFT(fv,fw,fx,fu)
         } else {
            if (u < x) a=u; else b=u;
            if (fu <= fw || w == x) {
               v=w;
               w=u;
               fv=fw;
               fw=fu;
            } else if (fu <= fv || v == x || v == w) {
               v=u;
               fv=fu;
            }
         }
      }
      fn_DisplayError("The maximum number of iterations in brent() is reached before convergence");
      *xmin=x;
      return fx;
   }
   #undef CGOLD
   #undef ZEPS
   #undef SHFT
   #undef SIGN

#else  //Default to CGII_OPTIMIZATION

   static int ncom;
   static double *pcom=NULL, *xicom=NULL, (*nrfunc)(double [], int); //nrfunc(), pcom, ncom, and xicom will be used by f1dim() and df1dim().
   static void (*nrdfun)(double [], double [], int, double (*func)(double [], int), double *, double);
   static void dlinmin(double p[], double xi[], int n, double *fret, double tol_dbrent, double itmax_dbrent, double *grdh_p, double (*func)(double [], int), void (*dfunc)(double [], double [], int, double (*func)(double [], int), double *, double)) {
      //Outputs:
      //  p[0, ..., n-1]:  a returned and reset value.
      //  xi[0, ..., n-1]:  a value repaced by the actual vector displacement that p was moved.
      //  fret:  the value of func at the returned location p.
      //Inputs:
      //  p[0, ..., n-1]:  a given point.
      //  xi[0, ..., n-1]:  a given multidimensional direction.
      //  n:  the dimension of p and xi.
      //  func():  the objective function.
      //  dfunc(): the gradient function computing the numerical gradient.  In the form of gradcd_gen() in cstz.c.

      int j;
      double xx,xmin,fx,fb,fa,bx,ax;

      ncom=n;
      pcom = tzMalloc(n, double);
      xicom = tzMalloc(n, double);
      nrfunc=func;
      nrdfun=dfunc;
      for (j=n-1;j>=0;j--) {
         pcom[j]=p[j];
         xicom[j]=xi[j];
      }
      ax=0.0;
      xx=1.0;
      mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
      *fret=dbrent(ax,xx,bx,f1dim, df1dim, grdh_p, tol_dbrent, itmax_dbrent, &xmin);
      for (j=n-1;j>=0;j--) {
         xi[j] *= xmin;
         p[j] += xi[j];
      }
      tzDestroy(xicom);
      tzDestroy(pcom);
   }


   //=== Used by dlinmin() only;
   #define ZEPS 1.0e-10
   #define MOV3(a,b,c, d,e,f)   {(a)=(d);(b)=(e);(c)=(f);}
   #define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
   static double dbrent(double ax, double bx, double cx, double (*f)(double), double (*df)(double, double *), double *grdh_p, double tol_dbrent, double itmax_dbrent, double *xmin) {
      int iter,ok1,ok2;
      double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
      double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

      a=(ax < cx ? ax : cx);
      b=(ax > cx ? ax : cx);
      x=w=v=bx;
      fw=fv=fx=(*f)(x);
      dw=dv=dx=(*df)(x, grdh_p);
      for (iter=1;iter<=itmax_dbrent;iter++) {
         xm=0.5*(a+b);
         tol1=tol_dbrent*fabs(x)+ZEPS;
         tol2=2.0*tol1;
         if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
            *xmin=x;
            return fx;
         }
         if (fabs(e) > tol1) {
            d1=2.0*(b-a);
            d2=d1;
            if (dw != dx) d1=(w-x)*dx/(dx-dw);
            if (dv != dx) d2=(v-x)*dx/(dx-dv);
            u1=x+d1;
            u2=x+d2;
            ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
            ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
            olde=e;
            e=d;
            if (ok1 || ok2) {
               if (ok1 && ok2)
                  d=(fabs(d1) < fabs(d2) ? d1 : d2);
               else if (ok1)
                  d=d1;
               else
                  d=d2;
               if (fabs(d) <= fabs(0.5*olde)) {
                  u=x+d;
                  if (u-a < tol2 || b-u < tol2)
                     d=SIGN(tol1,xm-x);
               } else {
                  d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
               }
            } else {
               d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
            }
         } else {
            d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
         }
         if (fabs(d) >= tol1) {
            u=x+d;
            fu=(*f)(u);
         } else {
            u=x+SIGN(tol1,d);
            fu=(*f)(u);
            if (fu > fx) {
               *xmin=x;
               return fx;
            }
         }
         du=(*df)(u, grdh_p);
         if (fu <= fx) {
            if (u >= x) a=x; else b=x;
            MOV3(v,fv,dv, w,fw,dw)
            MOV3(w,fw,dw, x,fx,dx)
            MOV3(x,fx,dx, u,fu,du)
         } else {
            if (u < x) a=u; else b=u;
            if (fu <= fw || w == x) {
               MOV3(v,fv,dv, w,fw,dw)
               MOV3(w,fw,dw, u,fu,du)
            } else if (fu < fv || v == x || v == w) {
               MOV3(v,fv,dv, u,fu,du)
            }
         }
      }
      fn_DisplayError("The maximum number of iterations in dbrent() is reached before convergence");
      return 0.0;
   }
   #undef ZEPS
   #undef MOV3
   #undef SIGN

   //=== Used by dlinmin() and dbrent() only;
   static double df1dim(double x, double *grdh_p) {
      int j;
      double df1=0.0;
      double *xt,*df;

      xt = tzMalloc(ncom, double);
      df = tzMalloc(ncom, double);
      for (j=ncom-1;j>=0;j--) xt[j]=pcom[j]+x*xicom[j];
      (*nrdfun)(df, xt, ncom, nrfunc, grdh_p, nrfunc(xt, ncom));
      //===================  WARNING ======================
      //We use 0.0 because the current gradient function gradcd_gen() in cstz.c do not use this function value.  A more
      //  sophisticated central gradient method would require this function value, and therefore we must pass
      //  nrfunc(xt, ncom) instead of 0.0.  TZ, September 2003.
      //===================  WARNING ======================
      for (j=ncom-1;j>=0;j--) df1 += df[j]*xicom[j];
      tzDestroy(df);
      tzDestroy(xt);
      return df1;
   }

#endif



static double f1dim(double x) {
   //Collapsing to one dimension line search, used by limin() or dlimin().
   int j;
   double f,*xt=NULL;

   xt = tzMalloc(ncom, double);
   for (j=ncom-1;j>=0;j--) xt[j]=pcom[j]+x*xicom[j];
   f=(*nrfunc)(xt, ncom);
   tzDestroy(xt);
   return f;
}


#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d)  {(a)=(b);(b)=(c);(c)=(d);}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (*func)(double)) {
   double ulim,u,r,q,fu,dum, tmpd;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
           (2.0*SIGN((tmpd=fabs(q-r))>TINY ? tmpd : TINY,q-r));   //Original: (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef SIGN





//-------------------
// My own functions.
//-------------------
//=== Computing Norm2 of dv.
static double ftd_norm2(double *vnew_p, double *vold_p, int _n) {
   int _i;
   double dtheta=0.0,  //Cumulative.
          tmpd;

   for (_i=_n-1; _i>=0; _i--) {
      tmpd = vnew_p[_i] - vold_p[_i];
      dtheta += square(tmpd);
   }

   return ( sqrt(dtheta) );
}

//=== Computing the inner product of x and y.
static double ftd_innerproduct(double *x, double *y, int _n) {
   int _i;
   double a = 0.0;   //Cumulative.
   for (_i=_n-1; _i>=0; _i--)  a += x[_i] * y[_i];    //a += (*x++) * (*y++);  Be aware that this alternative maybe too fancy.
   return (a);
}




//=== Extern function to be accessed by other C files.
void congradmin_SetPrintFile(char *filename) {
   if (!filename)   sprintf(filename_sp3vecs, "outdata5congradmin.prn");  //Default filename.
   else {
      strcpy(filename_sp3vecs, filename);
      //filename_sp3vecs[STRLEN-1] = '\0';  //The end of the string is set to NUL to prevent it from be a non-string.
   }
}



//void congradmin_SetPrintFile(FILE *fptr_sp) {
//   fptr_interesults = fptr_sp;
//}

//void congradmin_SetPrintFile_db(FILE *fptr_sp) {
//   fptr_interesults_db = fptr_sp;
//}


#undef STRLEN
