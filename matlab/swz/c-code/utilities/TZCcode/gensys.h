/*******************************************************************
 * [G1,C,impact,fmat,fwt,ywt,gev,eu]=gensys(g0,g1,c,psi,pi,div)
 *
 * System given as
 *         g0*y(t)=g1*y(t-1)+c+psi*z(t)+pi*eta(t),
 * with z an exogenous variable process and eta being endogenously determined
 * one-step-ahead expectational errors.  Returned system is
 *        y(t)=G1*y(t-1)+C+impact*z(t)+ywt*inv(I-fmat*inv(L))*fwt*z(t+1) .
 * If z(t) is i.i.d., the last term drops out.
 * If div or stake is omitted from argument list, a div>1 or stake>1 is calculated.
 * eu(1)=1 for existence, eu(2)=1 for uniqueness.  eu(1)=-1 for
 * existence only with not-serially correlated z(t); eu=[-2,-2] for coincident zeros.
 *
 * g0, g1:  n-by-n matrices.
 * c:  n-by-1 constant terms.
 * z(t):  m-by-1 vector of exogenous residuals where m < n.
 * psi:  n-by-m matrix.
 * eta(t):  h-by-1 vector of expectational (endogenous) errors.
 * pi:  n-by-h matrix.
 * div: a real number dividing stable and unstable roots..  If < 1.0, a div>1.0 is calculated mechanically.
 *-------
 * G1 or Theta_dm:  n-by-n matrices.
 * C:  n-by-1 vector of constant terms.
 * impact:  n-by-m matrix.
 * gev:  n-by-2 z vector of stacked generalized eigenvalues where gev(;,2) ./ gev(:,1) = eig(g0, g1).
 * ywt:  n-by-nunstab z matrix of possible complex numbers.  Initialized to NULL and dynamically allocated.
 * fmat: nunstab-by-nunstab z matrix where nunstab is the number of non-stable roots.
 * fwt:  nunstab-by-m z matrix.
********************************************************************/

#ifndef __GENSYS_H__
   #define __GENSYS_H__

   #include "tzmatlab.h"
   //#include "fn_filesetup.h"      //For DDDDebugging purpose.

   #define REALSMALL 1e-7
   //#define PRINTWARNINGofSUNSPOT

   typedef struct TSgensys_tag {
           //=== Output arguments.
           TSdmatrix *Theta_dm;  //n-by-n.
           TSdvector *c_dv;   //n-by-1.
           TSdmatrix *Impact_dm;  //n-by-m.
           TSdzmatrix *Fmat_dzm;   //nunstab-by-nunstab z matrix.  Initialized to NULL and will be dynamically allocated whenever gensys() is called.
           TSdzmatrix *Fwt_dzm;    //nunstab-by-m z matrix of possible complex numbers.  Initialized to NULL and dynamically allocated.
           TSdzmatrix *Ywt_dzm;    //n-by-nunstab z matrix of possible complex numbers.  Initialized to NULL and dynamically allocated.
           TSdzmatrix *Gev_dzm;  //n-by-2 z matrix of possible complex numbers.
           TSivector *eu_iv;   //2-by-1.
           //=== Function itself.
           int (*gensys)(struct TSgensys_tag *, void *);
           //=== Input arguments, which are all intialized to 0.0 and whose flags are set to M_GE.
           TSdmatrix *G0_dm;  //n-by-n.
           TSdmatrix *G1_dm;  //n-by-n.
           TSdvector *c0_dv;  //n-by-1.
           TSdmatrix *Psi_dm; //n-by-m.
           TSdmatrix *Pi_dm;  //n-by-k whtere k is the number of expectational errors.
           double div;  //Real number dividing stable and unstable roots..  If < 1.0, a div>1.0 is calculated mechanically.
   } TSgensys;
   //
   typedef int TFlinratexp(struct TSgensys_tag *, void *);  //For linear rational expectations models.

   struct TSgensys_tag *CreateTSgensys(TFlinratexp *func, const int _n, const int _m, const int _k, const double div);
   struct TSgensys_tag *DestroyTSgensys(struct TSgensys_tag *gensys_ps);
   int gensys_sims(struct TSgensys_tag *gensys_ps, void *dummy_ps);
#endif

