/*
* Copyright (C) 2008-2009 Dynare Team
*
* This file is part of Dynare.
*
* Dynare is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* Dynare is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
*/

/* derived from c++kalman_filter library by O. Kamenik */

/*****************
The file is divided into two parts, routines for multivariate filtering
and smoothing and univariate ones. The only reason for keeping these in one
file, is that the univariate routines use code bits from multivariate
routines.
*****************/

#include "kalman.h"
#include "ts_exception.h"

#include "cppblas.h"
//#include "cpplapack.h"

#include <math.h> 
#include <float.h> 
#include <cmath> 

/*****************
We delete everything which is not |NULL|. This is important, since
it may happen that the reults are not filled completely. (For instance
in the beginning or if error eccurred.)
*****************/

FilterResults::~FilterResults()
  {
  for(unsigned int i= 0;i<Finv.size();i++){
    if(Finv[i])
      delete Finv[i];
    if(v[i])
      delete v[i];
    if(L[i])
      delete L[i];
    if(a[i])
      delete a[i];
    if(P[i])
      delete P[i];
    }
  }

  /*****************
  We delete what is allocated and is in the way of the new data. Then
  set the new data as copies.
*****************/
void 
FilterResults::set(int t,const PLUFact&FFinv,const Vector&vv,
                   const GeneralMatrix&LL,const Vector&aa,
                   const GeneralMatrix&PP,double ll)
  {
  TS_RAISE_IF(t<1||t> (int)Finv.size()+1,
    "Wrong time for FilterResults::set");
  
  int tm= t-1;
  if(Finv[tm])
    delete Finv[tm];
  if(v[tm])
    delete v[tm];
  if(L[tm])
    delete L[tm];
  if(a[tm])
    delete a[tm];
  if(P[tm])
    delete P[tm];
  
  if(t> maxt)
    maxt= t;
  
  Finv[tm]= new PLUFact(FFinv);
  v[tm]= new Vector(vv);
  L[tm]= new GeneralMatrix(LL);
  a[tm]= new Vector(aa);
  P[tm]= new GeneralMatrix(PP);
  loglik[tm]= ll;
  }

const PLUFact&FilterResults::getFInverse(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getFInverse");
  return*(Finv[t-1]);
  }

const Vector&FilterResults::getV(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getV");
  return*(v[t-1]);
  }

const GeneralMatrix&FilterResults::getL(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getL");
  return*(L[t-1]);
  }

const Vector&FilterResults::getA(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getA");
  return*(a[t-1]);
  }

const GeneralMatrix&FilterResults::getP(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getP");
  return*(P[t-1]);
  }

  /*****************
  This adds all the log likelihoods for all periods. If some periods
  in the results have not been set, these are initialized to zeros and
  thus this method is pretty safe but may not be if the likelihood tends to
  be  far lower or higher than 0.
*****************/

double 
FilterResults::getLogLikelihood()const
  {
  double res= 0.0;
  for(unsigned int i= 0;i<loglik.size();i++)
    res+= loglik[i];
  return res;
  }

double 
FilterResults::getLogLikelihood(int start)const
  {
  double res= 0.0;
  for(unsigned int i= start;i<loglik.size();i++)
    res+= loglik[i];
  return res;
  }

double 
FilterResults::getLogLikelihood(std::vector<double>* vloglik)const
  {
  double res= 0.0;
  for(unsigned int i= 0;i<loglik.size();i++)
    res+= loglik[i];
  *vloglik= loglik;
  return res;
  }

double 
FilterResults::getLogLikelihood(int start,std::vector<double>* vloglik)const
  {
  double res= 0.0;
  for(unsigned int i= start;i<loglik.size();i++)
    res+= loglik[i];
  *vloglik= loglik;
  return res;
  }

DiffuseFilterResults::~DiffuseFilterResults()
  {
  for(unsigned int i= 0;i<L_1.size();i++){
    if(L_1[i])
      delete L_1[i];
    if(Pinf[i])
      delete Pinf[i];
    if(F_2[i])
      delete F_2[i];
    }
  }


bool 
DiffuseFilterResults::isFinfRegular(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for DiffuseFilterResults::isFinfRegular");
  return Finf_reg[t-1];
  }



bool 
DiffuseFilterResults::isPinfZero(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for DiffuseFilterResults::isPinfZero");
  return Pinf_zero[t-1];
  }

  /*****************
  Here we raise an error on attempt to retrieve inverse of $F_\infty$
  for a period in which it was not regular. Caller has to call
  |isFinfRegular| first.
*****************/

const PLUFact&DiffuseFilterResults::getFinfInverse(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for DiffuseFilterResults::getFinfInverse");
  TS_RAISE_IF(!isFinfRegular(t),
    "Finf not regular in the period in DiffuseFilterResults::getFinfInverse");
  return getFInverse(t);
  }

  /*****************
  Here we issue an error on attempt to retrieve inverse of $F_*$ for a
  period when $F_\infty$ was regular.
*****************/

const PLUFact&DiffuseFilterResults::getFstarInverse(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for DiffuseFilterResults::getFstarInverse");
  TS_RAISE_IF(isFinfRegular(t),
    "Finf not zero in the period in DiffuseFilterResults::getFstarInverse");
  return getFInverse(t);
  }

  /*****************
  This should be called only when $F_\infty$ was regular, we raise an
  error otherwise.
*****************/

const GeneralMatrix&DiffuseFilterResults::getF2(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getF2");
  TS_RAISE_IF(!isFinfRegular(t),
    "Finf not regular in the period in DiffuseFilterResults::getF2");
  return*(F_2[t-1]);
  }

  /*****************
  This should be called only when $F_\infty$ was regular, we raise an
  error otherwise.
*****************/

const GeneralMatrix&DiffuseFilterResults::getL1(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getL1");
  TS_RAISE_IF(!isFinfRegular(t),
    "Finf not regular in the period in DiffuseFilterResults::getL1");
  return*(L_1[t-1]);
  }

  /*****************
  The $P_\infty$ should be retrieved only if it is not zero, (these
  are all diffuse periods)
*****************/

const GeneralMatrix&DiffuseFilterResults::getPinf(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getPinf");
  TS_RAISE_IF(isPinfZero(t),
    "Pinf is zero in the period in DiffuseFilterResults::getPinf");
  return*(L_1[t-1]);
  }

  /*****************
  This sets the diffuse results for diffuse periods when
  $F_\infty$ is regular. Note that in this case the inverse
  $F_\infty^{-1}$ is stored as $F^{-1}$ of |FilterResults|, and thus is
  returned by call |getFinfInverse|. Also, $L^{(0)}$ is stored as $L$ of
  |FilterResults|, and retrieved by |getL0()| which is equivalent to
  |FilterResults::getL()| (see |@<|DiffuseFilterResults| class
  declaration@>|).
*****************/

void DiffuseFilterResults::set(int t,const PLUFact&FFinv,const GeneralMatrix&FF_2,
                               const Vector&vv,const GeneralMatrix&LL_0,
                               const GeneralMatrix&LL_1,const Vector&aa,
                               const GeneralMatrix&PPstar,const GeneralMatrix&PPinf,
                               double ll)
  {
  FilterResults::set(t,FFinv,vv,LL_0,aa,PPstar,ll);
  
  int tm= t-1;
  if(L_1[tm])
    delete L_1[tm];
  if(Pinf[tm])
    delete Pinf[tm];
  if(F_2[tm])
    delete F_2[tm];
  
  L_1[tm]= new GeneralMatrix(LL_1);
  Pinf[tm]= new GeneralMatrix(PPinf);
  F_2[tm]= new GeneralMatrix(FF_2);
  Finf_reg[tm]= true;
  Pinf_zero[tm]= false;
  }

  /*****************
  This sets the diffuse results for diffuse period when the $F_\infty$
  is zero. We do not set $F^{(2)}$, and we set $F_*^{-1}$ instead of
  $F_\infty^{-1}$.
*****************/

void DiffuseFilterResults::set(int t,const PLUFact&FFstarinv,const Vector&vv,
                               const GeneralMatrix&LL_0,const Vector&aa,
                               const GeneralMatrix&PPstar,const GeneralMatrix&PPinf,
                               double ll)
  {
  FilterResults::set(t,FFstarinv,vv,LL_0,aa,PPstar,ll);
  
  int tm= t-1;
  if(Pinf[tm])
    delete Pinf[tm];
  Pinf[tm]= new GeneralMatrix(PPinf);
  
  Finf_reg[tm]= false;
  Pinf_zero[tm]= false;
  }

  /*****************
  This returns number of initial diffuse periods (those having
  $P_\infty$ non-zero)
*****************/

int DiffuseFilterResults::getDiffusePeriods()const
  {
  int d= maxt;
  while(d> 1&&isPinfZero(d))
    d--;
  return d;
  }


;

SmootherResults::~SmootherResults()
  {
  for(unsigned int i= 0;i<alpha.size();i++){
    if(alpha[i])
      delete alpha[i];
    if(eta[i])
      delete eta[i];
    if(V[i])
      delete V[i];
    }
  }

;

void SmootherResults::set(int t,const Vector&aalpha,const Vector&eeta,
                          const GeneralMatrix&VV)
  {
  TS_RAISE_IF(t<1||t> (int)alpha.size()+1,
    "Wrong time for SmootherResults::set");
  if(t<mint)
    mint= t;
  int tm= t-1;
  if(alpha[tm])
    delete alpha[tm];
  if(eta[tm])
    delete eta[tm];
  if(V[tm])
    delete V[tm];
  
  alpha[tm]= new Vector(aalpha);
  eta[tm]= new Vector(eeta);
  V[tm]= new GeneralMatrix(VV);
  }

  /*****************
  This takes a |SmootherResults| coming from a univariate filter with
  a given number of observations at one time. This number of
  observations becomes a periodicity in the univariate
  |SmootherResults|. If this perioidicty is 10, we take data for
  $t=10,20,30,\ldots$.)
*****************/

void SmootherResults::import(const SmootherResults&sres,int period)
  {
  TS_RAISE_IF(period*alpha.size()!=sres.alpha.size(),
    "Results lengths not compatible with period in SmootherResults::import");
  TS_RAISE_IF(sres.mint!=1,
    "Results not finished in SmootherResults::import");
  for(unsigned int tm= 0;tm<alpha.size();tm++){
    if(alpha[tm])
      delete alpha[tm];
    if(eta[tm])
      delete eta[tm];
    if(V[tm])
      delete V[tm];
    alpha[tm]= new Vector((const Vector&)*sres.alpha[(tm+1)*period-1]);
    eta[tm]= new Vector((const Vector&)*sres.eta[(tm+1)*period-1]);
    V[tm]= new GeneralMatrix((const GeneralMatrix&)*sres.V[(tm+1)*period-1]);
    }
  
  mint= 1;
  }

  /*****************
  This saves |alpha| to the given matrix.
*****************/

void 
SmootherResults::exportAlpha(GeneralMatrix&out)const
  {
  TS_RAISE_IF(mint> 1,
    "Results not finished in SmootherResults::exportAlpha");
  TS_RAISE_IF(out.numCols()!=(int)alpha.size(),
    "Wrong number of columns in SmootherResults::exportAlpha");
  TS_RAISE_IF(alpha[0]->length()!=out.numRows(),
    "Wrong number of rows in SmootherResults::exportAlpha");
  for(unsigned int tm= 0;tm<alpha.size();tm++)
    {
    Vector outi(out,tm);
    outi= (const Vector&)(*alpha[tm]);
    }
  }

  /*****************
  This saves |eta| to the given matrix.
*****************/

void 
SmootherResults::exportEta(GeneralMatrix&out)const
  {
  TS_RAISE_IF(mint> 1,
    "Results not finished in SmootherResults::exportEta");
  TS_RAISE_IF(out.numCols()!=(int)eta.size(),
    "Wrong number of columns in SmootherResults::exportEta");
  TS_RAISE_IF(eta[0]->length()!=out.numRows(),
    "Wrong number of rows in SmootherResults::exportEta");
  for(unsigned int tm= 0;tm<eta.size();tm++)
    {
    Vector outi(out,tm);
    outi= (const Vector&)(*eta[tm]);
    }
  }

  /*****************
  This saves $V$ to the given two dimensional matrix. We store the $V$
  matrices one by one columnwise. The storage corresponds to Matlab
  storage of three dimensional matrices.
*****************/

void SmootherResults::exportV(GeneralMatrix&out)const
  {
  TS_RAISE_IF(mint> 1,
    "Results not finished in SmootherResults::exportV");
  int m= V[0]->numRows();
  TS_RAISE_IF(out.numCols()!=(int)V.size()*m,
    "Wrong number of columns in SmootherResults::exportV");
  TS_RAISE_IF(m!=out.numRows(),
    "Wrong number of rows in SmootherResults::exportV");
  for(unsigned int tm= 0;tm<V.size();tm++)
    {
    GeneralMatrix outi(out,0,tm*m,m,m);
    outi= (const GeneralMatrix&)(*V[tm]);
    }
  }

BasicKalmanTask::BasicKalmanTask(const GeneralMatrix&d,const GeneralMatrix&ZZ,
                                 const GeneralMatrix&HH,const GeneralMatrix&TT,
                                 const GeneralMatrix&RR,const GeneralMatrix&QQ,
                                 const StateInit&init_state)
                                 :  // ssf(Z,H,T,R,Q),
data(d), Zt(*(new ConstGeneralMatrix(ZZ))), 
Ht(*(new ConstGeneralMatrix(HH))), 
Tt(*(new ConstGeneralMatrix(TT))), 
Rt(*(new ConstGeneralMatrix(RR))), 
Qt(*(new ConstGeneralMatrix(QQ))),
init(init_state)
  {
  TS_RAISE_IF(d.numRows()!=Zt.numRows(),
    "Data not compatible with BasicKalmanTask constructor");
  //  TS_RAISE_IF(ssf.m!=init.getM(),
  //    "State initialization not compatible with SSF in KalmanTask constructor");
  }

BasicKalmanTask::BasicKalmanTask(const GeneralMatrix&d,const ConstGeneralMatrix&ZZ,
                                 const ConstGeneralMatrix&HH,const ConstGeneralMatrix&TT,
                                 const ConstGeneralMatrix&RR,const ConstGeneralMatrix&QQ,
                                 const StateInit&init_state)
                                 :  // ssf(Z,H,T,R,Q),
data(d), Zt(ZZ), Ht(HH), Tt(TT), Rt(RR), Qt(QQ),init(init_state)
  {
  TS_RAISE_IF(d.numRows()!=Zt.numRows(),
    "Data not compatible with BasicKalmanTask constructor");
  //  TS_RAISE_IF(ssf.m!=init.getM(),
  //    "State initialization not compatible with SSF in KalmanTask constructor");
  }


BasicKalmanTask::~BasicKalmanTask()
  {
  if (&Zt);
    delete &Zt;
  if (&Ht);
    delete &Ht;
  if (&Tt);
    delete &Tt;
  if (&Rt);
    delete &Rt;
  if (&Qt);
    delete &Qt;
  }

KalmanTask::KalmanTask(const GeneralMatrix&d,const GeneralMatrix&Z,
                       const GeneralMatrix&H,const GeneralMatrix&T,
                       const GeneralMatrix&R,const GeneralMatrix&Q,
                       const StateInit&init_state)
                       :ssf(Z,H,T,R,Q),
                       data(d),
                       init(init_state)
  {
  TS_RAISE_IF(d.numRows()!=Z.numRows(),
    "Data not compatible with SSF in KalmanTask constructor");
  TS_RAISE_IF(ssf.m!=init.getM(),
    "State initialization not compatible with SSF in KalmanTask constructor");
  }

KalmanTask::KalmanTask(const GeneralMatrix&d,const TMatrix&Z,
                       const TMatrix&H,const TMatrix&T,
                       const TMatrix&R,const TMatrix&Q,
                       const StateInit&init_state)
                       :ssf(Z,H,T,R,Q),
                       data(d),
                       init(init_state)
  {
  TS_RAISE_IF(d.numRows()!=Z.numRows(),
    "Data not compatible with SSF in KalmanTask constructor");
  TS_RAISE_IF(ssf.m!=init.getM(),
    "State initialization not compatible with SSF in KalmanTask constructor");
  }

  /*****************
  This is a public interface to mechanism that filters data and returns 
  the data loglikelihood. In addition,
  it returns the number of periods |per| successfully filtered. If the
  filter failed to filter entire data, than the resulting |per| is less
  than a number of columns in the |data|. This might be caused by
  singular $F$ matrix, or by singular and non-zero $F_\infty$
  matrix. Further it returns a number od diffuse initial periods, in
  case of non-diffuse initialization, zero is returned.
*****************/
double 
KalmanTask::filter(int&per,int&d)const
  {
  if(!init.isDiffuse())
    {
    FilterResults fres(data.numCols());
    filterNonDiffuse(init.getA(),init.getPstar(),1,fres);
    d= 0;
    per= fres.getMaxT();
    if(fres.hasFinished())
      return fres.getLogLikelihood();
    }
  else
    {
    DiffuseFilterResults fres(data.numCols());
    filterDiffuse(init.getA(),init.getPstar(),init.getPinf(),
      1,fres);
    d= fres.getDiffusePeriods();
    per= fres.getMaxT();
    if(fres.hasFinished())
      return fres.getLogLikelihood();
    }
  return 0.0;
  }

double 
KalmanTask::filter(int&per,int&d, int start, std::vector<double>* vll)const
  {
  if(!init.isDiffuse())
    {
    FilterResults fres(data.numCols());
    filterNonDiffuse(init.getA(),init.getPstar(),1,fres);
    d= 0;
    per= fres.getMaxT();
    if(fres.hasFinished())
      return fres.getLogLikelihood(start, vll);
    }
  else
    {
    DiffuseFilterResults fres(data.numCols());
    filterDiffuse(init.getA(),init.getPstar(),init.getPinf(),
      1,fres);
    d= fres.getDiffusePeriods();
    per= fres.getMaxT();
    if(fres.hasFinished())
      return fres.getLogLikelihood(start, vll);
    }
  return 0.0;
  }

double 
BasicKalmanTask::filter(int&per,int&d, int start, std::vector<double>* vll)const
  {
  d= 0;
  per= vll->size() ;
  return filterNonDiffuse(init.getA(),init.getPstar(), start, vll);
  }

/*****************
  This is public interface that runs a filter followed by a smoother. 
  In addition to things returned by |KalmanTask::filter|, it fills also 
  |SmootherResults|, which must be initialized to the number of columns 
  of the |data|.
*****************/
double 
KalmanTask::filter_and_smooth(SmootherResults&sres,
                              int&per,int&d)const
  {
  if(!init.isDiffuse())
    {
    FilterResults fres(data.numCols());
    filterNonDiffuse(init.getA(),init.getPstar(),1,fres);
    d= 0;
    per= fres.getMaxT();
    if(fres.hasFinished())
      {
      smootherNonDiffuse(fres,sres);
      return fres.getLogLikelihood();
      }
    }
  else
    {
    DiffuseFilterResults fres(data.numCols());
    filterDiffuse(init.getA(),init.getPstar(),init.getPinf(),
      1,fres);
    d= fres.getDiffusePeriods();
    per= fres.getMaxT();
    if(fres.hasFinished())
      {
      smootherDiffuse(fres,sres);
      return fres.getLogLikelihood();
      }
    }
  return 0.0;
  }


/*****************
  This runs a Basic non-diffuse filter with the given $t$, $a_t$ and
  $P_t$. It fills the |FilterResults|.
  
    First we check that the passed $P_t$ is positive definite by checking
    that it has not negative diagonal and is symmetric diagonally
    dominant. This is not equivalent to positive definitness but it excludes
    ``much'' of indefinite matrices. This check is important since this
    routine is called from a diffuse filter and it is possible due to a
    wrong guess/calculation of a number of diffuse periods that numerical
    instability and roundoff errors make the matrix $P_*$ broken.
    
      Then we cycle until the end of period and perform a classical Kalman
      filter operations.
*****************/
double 
BasicKalmanTask::filterNonDiffuse(const Vector&a,const GeneralMatrix&P,
                                  int start, std::vector<double>* vll)const
  {
  double loglik=0;
  Vector at(a);
  GeneralMatrix Pt(P);
//  GeneralMatrix PtZeros(Pt.numRows(), Pt.numCols());
//  PtZeros.zeros();
  if(TSUtils::hasNegativeDiagonal(Pt)||!TSUtils::isSymDiagDominant(Pt))
    return 0.0;
  const int m=Pt.numRows();
  const int n=Zt.numRows();
  int inc =1;
  const int rcols= Rt.numCols();
  GeneralMatrix Ft (Ht.numRows(), Ht.numCols()  );
  GeneralMatrix Lt(Tt);
  GeneralMatrix PtLttrans(m,m);
  GeneralMatrix Mt(m,n);
  GeneralMatrix Kt(m,n);
  GeneralMatrix QtRttrans(rcols,Rt.numRows());
//  PLUFact Ftinv(Ft);

  bool isTunit=0;// Tt->isUnit();
  bool isQzero= Qt.isZero();
  bool isRzero= Rt.isZero();
  const double alpha=1.0;
  const double neg_alpha=-1.0;
  const double omega=0.0;
  Vector vt(n);
  Vector atsave(m);//((const Vector&)at);
  int t= 1; 
  int nonsteady=1;
  for(;t<=data.numCols()&&nonsteady;++t)
    {
    ConstVector yt(data,t-1);
    
    /*****************
    This calculates  $$v_t = y_t - Z_t*a_t.$$
    *****************/   
//    Vector vt(yt);
    vt=yt;
//    Zt.multsVec(vt,at);
	  BLAS_dgemv("N",  &n,  &m, &neg_alpha, Zt.base(), &n, at.base(), 
        &inc, &alpha, vt.base(), &inc);
    
    /*****************
    This calculates $$F_t = Z_tP_tZ_t^T+H_t.$$
    *****************/   
//    GeneralMatrix Mt(Pt,Zt,"trans");
		BLAS_dgemm("N", "T", &m, &n, &m, &alpha, Pt.base(), &m,
				   Zt.base(), &n, &omega, Mt.base(), &m); 

//    GeneralMatrix Ft(Ht);
    Ft=Ht;
//    Ft.multAndAdd(Zt,ConstGeneralMatrix(Mt));
    // DGEMM: C := alpha*op( A )*op( B ) + beta*C,
		BLAS_dgemm("N", "N", &n, &n, &m, &alpha, Zt.base(), &n,
				   Mt.base(), &m, &alpha, Ft.base(), &n); 
    
    PLUFact Ftinv(Ft); //    Ftinv=Ft;
    if(!Ftinv.isRegular())
      return 0.0;
    
    /*****************
    This calculates $$K_t = T_tP_tZ_t^TF_t^{-1}.$$
    *****************/   
//    GeneralMatrix Kt(Tt,Mt);
		BLAS_dgemm("N", "N", &m, &n, &m, &alpha, Tt.base(), &m,
				   Mt.base(), &m, &omega, Kt.base(), &m); 
    Ftinv.multInvRight(Kt);
    
    /*****************
    This calculates  $$L_t = T_t-K_tZ_t.$$
    *****************/   
    //GeneralMatrix Lt(Tt);
    Lt=Tt;
    //Lt.multAndAdd(ConstGeneralMatrix(Kt),Zt,-1.0);
    // DGEMM: C := alpha*op( A )*op( B ) + beta*C,
		BLAS_dgemm("N", "N", &m, &m, &n, &neg_alpha, Kt.base(), &m,
				   Zt.base(), &n, &alpha, Lt.base(), &m); 
    
    
    /*****************
    Here we calc likelihood and store results.
    *****************/   
    double ll= calcStepLogLik(Ftinv,vt);
    //    fres.set(t,Ftinv,vt,Lt,at,Pt,ll);
    (*vll)[t-1]=ll;
    if (t>start) loglik+=ll;
    
    if(t<data.numCols())
      {
      /*****************
      This calculates $$a_{t+1} = T_ta_t + K_tv_t.$$
      *****************/   
      if(!isTunit)
        {
//        Vector atsave((const Vector&)at);
        atsave=at;
        Tt.multVec(0.0,at,1.0,atsave);
        }
      Kt.multVec(1.0,at,1.0,ConstVector(vt));
      
      /*****************
      This calculates $$P_{t+1} = T_tP_tL_t^T + R_tQ_tR_t^T.$$
      *****************/   
//    GeneralMatrix PtLttrans(Pt,Lt,"trans");
   // DGEMM: C := alpha*op( A )*op( B ) + beta*C,
  		BLAS_dgemm("N", "T", &m, &m, &m, &alpha, Pt.base(), &m,
				   Lt.base(), &m, &omega, PtLttrans.base(), &m); 
      if(!isTunit)
        {
//        Pt.zeros();
//        Pt=PtZeros;
//        Pt.multAndAdd(Tt,ConstGeneralMatrix(PtLttrans));
        // DGEMM: C := alpha*op( A )*op( B ) + beta*C,
		    BLAS_dgemm("N", "N", &m, &m, &m, &alpha, Tt.base(), &m,
				       PtLttrans.base(), &m, &omega, Pt.base(), &m); 
        }
      else
        {
        Pt= (const GeneralMatrix&)PtLttrans;
        }
      if(!isRzero&&!isQzero)
        {
//      GeneralMatrix QtRttrans(Qt,Rt,"trans");
        BLAS_dgemm("N", "T", &rcols, &m, &rcols, &alpha, Qt.base(), &rcols,
				   Rt.base(), &m, &omega, QtRttrans.base(), &rcols); 

 //       Pt.multAndAdd(Rt,ConstGeneralMatrix(QtRttrans));
        // DGEMM: C := alpha*op( A )*op( B ) + beta*C,
		    BLAS_dgemm("N", "N", &m, &m,  &rcols, &alpha, Rt.base(), &m,
				       QtRttrans.base(), &rcols, &alpha, Pt.base(), &m); 
        }
      
      }
    }
//  for(;t<=data.numCols();t++)
//    {
//    ConstVector yt(data,t-1);
//    }
  return loglik;
  }
  
  
  double 
    BasicKalmanTask::calcStepLogLik(const PLUFact&Finv,const Vector&v)
    {
    int p= Finv.numRows();
    Vector Finvv(v);
    Finv.multInvLeft(Finvv);
    double vFinvv= v.dot(Finvv);
    return-0.5*(p*log(2*M_PI)+Finv.getLogDeterminant()+vFinvv);
    }
  
  
    /*****************
    This runs a non-diffuse filter with the given $t$, $a_t$ and
    $P_t$. It fills the |FilterResults|.
    
      First we check that the passed $P_t$ is positive definite by checking
      that it has not negative diagonal and is symmetric diagonally
      dominant. This is not equivalent to positive definitness but it excludes
      ``much'' of indefinite matrices. This check is important since this
      routine is called from a diffuse filter and it is possible due to a
      wrong guess/calculation of a number of diffuse periods that numerical
      instability and roundoff errors make the matrix $P_*$ broken.
      
        Then we cycle until the end of period and perform a classical Kalman
        filter operations.
  *****************/
  void 
    KalmanTask::filterNonDiffuse(const Vector&a,const GeneralMatrix&P,
    int first,FilterResults&fres)const
    {
    Vector at(a);
    GeneralMatrix Pt(P);
    if(TSUtils::hasNegativeDiagonal(Pt)||!TSUtils::isSymDiagDominant(Pt))
      return;
    
    for(int t= first;t<=data.numCols();t++)
      {
      ConstVector yt(data,t-1);
      ConstGeneralMatrix Zt(((const TMatrix&)*ssf.Z)[t]);
      ConstGeneralMatrix Ht(((const TMatrix&)*ssf.H)[t]);
      ConstGeneralMatrix Tt(((const TMatrix&)*ssf.T)[t]);
      ConstGeneralMatrix Qt(((const TMatrix&)*ssf.Q)[t]);
      ConstGeneralMatrix Rt(((const TMatrix&)*ssf.R)[t]);
      bool isTunit= ssf.T->isUnit(t);
      bool isQzero= ssf.Q->isZero(t);
      bool isRzero= ssf.R->isZero(t);
      
      /*****************
      This calculates  $$v_t = y_t - Z_t*a_t.$$
      *****************/   
      Vector vt(yt);
      Zt.multsVec(vt,at);
      
      /*****************
      This calculates $$F_t = Z_tP_tZ_t^T+H_t.$$
      *****************/   
      GeneralMatrix Mt(Pt,Zt,"trans");
      GeneralMatrix Ft(Ht);
      Ft.multAndAdd(Zt,ConstGeneralMatrix(Mt));
      
      
      PLUFact Ftinv(Ft);
      if(!Ftinv.isRegular())
        return;
      
        /*****************
        This calculates $$K_t = T_tP_tZ_t^TF_t^{-1}.$$
      *****************/   
      GeneralMatrix Kt(Tt,Mt);
      Ftinv.multInvRight(Kt);
      
      /*****************
      This calculates  $$L_t = T_t-K_tZ_t.$$
      *****************/   
      GeneralMatrix Lt(Tt);
      Lt.multAndAdd(ConstGeneralMatrix(Kt),Zt,-1.0);
      
      
      /*****************
      Here we calc likelihood and store results.
      *****************/   
      double ll= calcStepLogLik(Ftinv,vt);
      fres.set(t,Ftinv,vt,Lt,at,Pt,ll);
      
      if(t<data.numCols())
        {
        /*****************
        This calculates $$a_{t+1} = T_ta_t + K_tv_t.$$
        *****************/   
        if(!isTunit)
          {
          Vector atsave((const Vector&)at);
          Tt.multVec(0.0,at,1.0,atsave);
          }
        Kt.multVec(1.0,at,1.0,ConstVector(vt));
        
        /*****************
        This calculates $$P_{t+1} = T_tP_tL_t^T + R_tQ_tR_t^T.$$
        *****************/   
        GeneralMatrix PtLttrans(Pt,Lt,"trans");
        if(!isTunit)
          {
          Pt.zeros();
          Pt.multAndAdd(Tt,ConstGeneralMatrix(PtLttrans));
          }
        else
          {
          Pt= (const GeneralMatrix&)PtLttrans;
          }
        if(!isRzero&&!isQzero)
          {
          GeneralMatrix QtRttrans(Qt,Rt,"trans");
          Pt.multAndAdd(Rt,ConstGeneralMatrix(QtRttrans));
          }
        }
      }
    }
  
/*****************
This runs a diffuse filter. Similarly as for
|KalmanTask::filterNonDiffuse| the filter is started with a given $t$,
$a_t$, $P_*$, and $P_\infty$ and stores the results to
|DiffuseFilterResults| |fres|.
    
This executes the diffuse multivariate filter period by period and if
the variance of states $P=P_*+\kappa P_\infty$ is finite for
$kappa->oo $, then we switch to |KalmanTask::filterNonDiffuse|. 

The switching has two reasons: 
The first is that the non-diffuse filter is computationally more efficient
(since it avoids multiplications of zero matrices). The second reason
is much more important. As $P_\infty$ approaches to zero, then
$F_\infty=Z P_\infty Z^T$ approaches to zero and might contain severe
roundoffs. All the operations employing its inverse, $F_\infty^{-1}$,
will commit very bad roundoff errors, and the results will become
unusable. That is why it is important to not only switch to
non-diffuse filter, but also to switch at the right period.

In theory, the period $d$ of switching is equal to a number of
(univariate) observations for which $F_\infty$ is regular. This is
because the regular $F_\infty=ZP_\infty Z^T$ conveys some information
to $P=P_*+\kappa P_\infty$. However, it is only a theoretical result;
in real floating point world it is difficult to recognize a regular
matrix in this process. Moreover, the $F_\infty$ might be singular and
still convey some information for the diffuse elements since it might
have non-zero rank.

In this implementation, we use the above idea with the following test
for regularity of $F_\infty$. $F_\infty$ is considered to be regular,
if its PLU factorization yields a condition number estimate less than
$10^10$. During the process it might happen that $P_\infty$ is
indefinite. In this case we correct it by setting its negative
eigenvalues to zero. So $F_\infty=ZP_\infty Z$ is always positive
semidefinite, so no tests for a sign of its determinant are
necessary. Further, the test for $F_\infty=0$ here is equivalent to an
exact match. This can be done since the roundoff errors are believed
to be eliminated during correcting the $P_\infty$ matrix, where not
only negative eigenvalues but also very small positive eigenvalues are
corrected to zeros. In neither case, this is if $F_\infty$ is regular
and still is non-zero, we raise end the filter. This error can be
recognized by |FilterResults::per| less than a number of periods.

This is just one of many ways, how to implement this non-continuous
algorithm. It is theoretically continuous (since the non-diffuse
periods having $P_\infty$ zero are covered by the branch where
$F_\infty=0$). However, it is computationally discontinuous, since the
calcs depend on when we switch to non-diffuse filter. Because of the
roundoff errors we are uncertain about the switch. An experience shows
that if we switch late, the results can be very bad due to roundoff
errors implied by late switch, if we switch too early, the results
might be wrong since we neglect some uncertainity.

Main decision point is |ndiff|. Whenever |ndiff<=0|, we consider
$P_\infty$ as zero and carry on as in non-diffuse filter.
*****************/   

void 
KalmanTask::filterDiffuse(const Vector&a,const GeneralMatrix&Pstar,
const GeneralMatrix&Pinf,int first,
DiffuseFilterResults&fres)const
  {
  Vector at(a);
  GeneralMatrix Ptstar(Pstar);
  GeneralMatrix Ptinf(Pinf);
  int ndiff= init.getNDiff();
  for(int t= first;t<=data.numCols();t++)
    {
    
    /*****************
    If $P_\infty$ is exactly zero, then we run the non-diffuse
    filter. The $P_\infty$ might become exactly zero by negative or zero
    |ndiff|, or by $P\infty$ definitness correction.    
    *****************/   
    if(TSUtils::isZero(Ptinf))
      {
      filterNonDiffuse(at,Ptstar,t,fres);
      return;
      }
    
    ConstVector yt(data,t-1);
    ConstGeneralMatrix Zt(((const TMatrix&)*ssf.Z)[t]);
    ConstGeneralMatrix Ht(((const TMatrix&)*ssf.H)[t]);
    ConstGeneralMatrix Tt(((const TMatrix&)*ssf.T)[t]);
    ConstGeneralMatrix Qt(((const TMatrix&)*ssf.Q)[t]);
    ConstGeneralMatrix Rt(((const TMatrix&)*ssf.R)[t]);
    bool isTunit= ssf.T->isUnit(t);
    bool isQzero= ssf.Q->isZero(t);
    bool isRzero= ssf.R->isZero(t);
    
    /*****************
    This calculates  $$v_t = y_t - Z_t*a_t.$$
    *****************/   
    Vector vt(yt);
    Zt.multsVec(vt,at);
    
    
    /*****************
    This calculates $$M_{*,t} = P_{*,t}Z_t^T.$$
    *****************/   
    GeneralMatrix Mtstar(Ptstar,Zt,"trans");
    
    /*****************
    This calculates $$F_{*,t} = Z_tP_{*,t}^T+H_t.$$
    *****************/   
    GeneralMatrix Ftstar(Ht);
    Ftstar.multAndAdd(Zt,ConstGeneralMatrix(Mtstar));
    
    /*****************
    This calculates $$M_{\infty,t} = P_{\infty,t}Z_t^T.$$
    *****************/   
    GeneralMatrix Mtinf(Ptinf,Zt,"trans");
    
    /*****************
    This calculates $$F_{\infty,t} = Z_tP_{\infty,t}Z_t^T.$$
    *****************/   
    GeneralMatrix Ftinf(Zt,ConstGeneralMatrix(Mtinf));
    
    
    PLUFact Ftinfinv(Ftinf);
    if(Ftinfinv.isRegular()&&Ftinfinv.getRcond()> 1.e-10)
      {
      ndiff-= ssf.p;
      
      /*****************
      We calculate all other matrices, and if we have not come to the end,
      also $a_{t+1}$, $P_{*,t+1}$ and $P_{\infty,t+1}$. If |ndiff<=0|, we
      set $P_{\infty,t+1}=0$. The matrix can be set to zero even if it is
      not positive semidefinite in the code correcting definitness of
      $P_\infty$.
      *****************/   
      
      /*****************
      This calculates $$F_t^{(2)} = -F_{\infty,t}^{-1}F_{*,t}F_{\infty,t}^{-1}.$$
      *****************/   
      GeneralMatrix Ft_2(Ftstar);
      Ftinfinv.multInvRight(Ft_2);
      Ftinfinv.multInvLeft(Ft_2);
      Ft_2.mult(-1.0);
      
      /*****************
      This calculates $$K_t^{(0)} = T_tM_{\infty,t}F_t^{(1)}.$$
      *****************/   
      GeneralMatrix Kt_0(Tt,Mtinf);
      Ftinfinv.multInvRight(Kt_0);
      
      /*****************
      This calculates $$K_t^{(1)} = T_t(M_{\infty,t}F_t^{(2)}+M_{*,t}F_t^{(1)}).$$
      *****************/   
      GeneralMatrix Kt_1(Mtstar);
      Ftinfinv.multInvRight(Kt_1);
      Kt_1.multAndAdd(Mtinf,Ft_2);
      if(!isTunit)
        Kt_1.multLeft(Tt);
      
      /*****************
      This calculates $$L_t^{(0)} = T_t-K_t^{(0)}Z_t.$$
      *****************/   
      GeneralMatrix Lt_0(Tt);
      Lt_0.multAndAdd(ConstGeneralMatrix(Kt_0),Zt,-1.0);
      
      /*****************
      This calculates $$L_t^{(1)} = -K_t^{(1)}Z_t.$$
      *****************/   
      GeneralMatrix Lt_1(Kt_1,Zt);
      Lt_1.mult(-1.0);
      
      /*****************
      This calculates log likelihood and store results
      *****************/   
      double ll= -0.5*(ssf.p*log(2*M_PI)+Ftinfinv.getLogDeterminant());
      fres.set(t,Ftinfinv,Ft_2,vt,Lt_0,Lt_1,at,Pstar,Pinf,ll);
      
      if(t<data.numCols())
        {
        if(!isTunit)
          {
          Vector atsave((const Vector&)at);
          Tt.multVec(0.0,at,1.0,atsave);
          }
        Kt_0.multVec(1.0,at,1.0,vt);
        
        /*****************
        This calculates
        $$P_{*,t+1} = T_t(P_{*,t}L_t^{(0)T}+P_{\infty,t}L_t^{(1)T})+R_tQ_tR_t^T.$$
        *****************/
        GeneralMatrix tmp(Ptstar,Lt_0,"trans");
        tmp.multAndAdd(Ptinf,Lt_1,"trans");
        if(!isTunit)
          Ptstar.mult(Tt,ConstGeneralMatrix(tmp));
        else
          Ptstar= (const GeneralMatrix&)tmp;
        if(!isQzero&&!isRzero)
          {
          GeneralMatrix QtRttrans(Qt,Rt,"trans");
          Ptstar.multAndAdd(Rt,ConstGeneralMatrix(QtRttrans));
          }
        
        /*****************
        We call |TSUtils::correctDefinitness| only if it has a negative
        diagonal or it is not diagonall dominant. We could call the routine in
        all any case, but it is costly.        
        *****************/
        if(TSUtils::hasNegativeDiagonal(Ptstar)||!TSUtils::isSymDiagDominant(Ptstar))
          TSUtils::correctDefinitness(Ptstar);
        
        /*****************
        This calculates $$P_{\infty,t+1} = T_tP_{\infty,t}L_t^{(0)T}.$$ Due
        to possible roundoff errors, the resulting matrix might not be
        symmetric, so we amend it by putting it to ${1\over
        2}(P_{\infty,t+1}+P_{\infty,t+1}^T)$.        
        *****************/
        if(!isTunit)
          Ptinf.multLeft(Tt);
        Ptinf.multRightTrans(Lt_0);
        TSUtils::correctSymmetricity(Ptinf);
        
        /*****************
        We check the semidefinitness of new $P_{\infty,t+1}$. If it is not,
        then the roundoff error is guilty for the mess and we have to correct
        the matrix to be semidefinite.        
        *****************/
        if(TSUtils::hasNegativeDiagonal(Ptinf)||!TSUtils::isSymDiagDominant(Ptinf))
          TSUtils::correctDefinitness(Ptinf);
        
        if(ndiff<=0)
          Ptinf.zeros();
        }
      }
    else if(TSUtils::isZero(Ftinf))
      {
      /*****************
      If $F_{*,t}$ is not regular, we return and the filter has not
      finished. The regularity is checked exactly.        
      *****************/
      PLUFact Ftstarinv(Ftstar);
      if(!Ftstarinv.isRegular())
        {
        return;
        }
      
      /*****************
      This calculates $$K_t^{(0)} = T_tM_{*,t}F_{*,t}^{-1}.$$
      *****************/   
      GeneralMatrix Kt_0(Tt,Mtstar);
      Ftstarinv.multInvRight(Kt_0);
      
      /*****************
      This calculates $$L_t^{(0)} = T_t-K_t^{(0)}Z_t.$$
      *****************/   
      GeneralMatrix Lt_0(Tt);
      Lt_0.multAndAdd(ConstGeneralMatrix(Kt_0),Zt,-1.0);
      
      
      /*****************
      This calculates log likelihood and store results
      *****************/   
      double ll= calcStepLogLik(Ftstarinv,vt);
      fres.set(t,Ftstarinv,vt,Lt_0,at,Ptstar,Ptinf,ll);
      
      
      if(t<data.numCols())
        {
        /*****************
        This calculates $$a_{t+1} = T_ta_t+K_t^{(0)}v_t.$$
        *****************/   
        if(!isTunit)
          {
          Vector atsave((const Vector&)at);
          Tt.multVec(0.0,at,1.0,atsave);
          }
        Kt_0.multVec(1.0,at,1.0,vt);
        
        /*****************
        This calculates $$P_{\infty,t+1} = T_tP_{\infty,t}T_t^T.$$
        *****************/   
        if(!isTunit)
          {
          GeneralMatrix PtinfTttrans(Ptinf,Tt,"trans");
          Ptinf.mult(Tt,ConstGeneralMatrix(PtinfTttrans));
          }
        if(TSUtils::hasNegativeDiagonal(Ptinf)||!TSUtils::isSymDiagDominant(Ptinf))
          TSUtils::correctDefinitness(Ptinf);
        
        /*****************
        This calculates $$P_{*,t+1} = T_tP_{*,t}L_t^{(0)T}+R_tQ_tR_t^T.$$   
        *****************/   
        GeneralMatrix PtstarLt_0trans(Ptstar,Lt_0,"trans");
        if(!isTunit)
          Ptstar.mult(Tt,ConstGeneralMatrix(PtstarLt_0trans));
        else
          Ptstar= (const GeneralMatrix&)PtstarLt_0trans;
        if(!isQzero&&!isRzero)
          {
          GeneralMatrix QtRttrans(Qt,Rt,"trans");
          Ptstar.multAndAdd(Rt,ConstGeneralMatrix(QtRttrans));
          }
        
        if(TSUtils::hasNegativeDiagonal(Ptstar)||!TSUtils::isSymDiagDominant(Ptstar))
          TSUtils::correctDefinitness(Ptstar);
        }
      }
    else
      {
      return;
      }
    }
  }
  
  
  /*****************
  This executes only one step of smoother non-diffuse step. It takes
  $r_t$, and $N_t$ and outputs $\alpha_t$, $V_t$ and $\eta_t$. The code is clear.
  *****************/
  void 
  KalmanTask::smootherNonDiffuseStep(int t,const FilterResults&fres,
    Vector&rt,GeneralMatrix&Nt,
    Vector&alphat,GeneralMatrix&Vt,
    Vector&etat)const
    {
    const PLUFact&Ftinv= fres.getFInverse(t);
    const Vector&vt= fres.getV(t);
    const GeneralMatrix&Lt= fres.getL(t);
    const Vector&at= fres.getA(t);
    const GeneralMatrix&Pt= fres.getP(t);
    ConstGeneralMatrix Zt(((const TMatrix&)*ssf.Z)[t]);
    ConstGeneralMatrix Qt(((const TMatrix&)*ssf.Q)[t]);
    ConstGeneralMatrix Rt(((const TMatrix&)*ssf.R)[t]);
    bool isQzero= ssf.Q->isZero(t);
    bool isRzero= ssf.R->isZero(t);
    
    
    /*****************
    Calculate $$\eta_t = Q_tR_t^Tr_t.$$
    *****************/
    etat.zeros();
    if(!isQzero&&!isRzero){
      Rt.multVecTrans(0.0,etat,1.0,rt);
      Vector etatsav((const Vector&)etat);
      Qt.multVec(0.0,etat,1.0,etatsav);
      }
    
    
    
    /*****************
    This calculates $$r_{t-1} = Z^T_tF_t^{-1}v_t + L^T_tr_t.$$
    *****************/
    Vector rtsav((const Vector&)rt);
    Lt.multVecTrans(0.0,rt,1.0,rtsav);
    Vector Ftinvvt(vt);
    Ftinv.multInvLeft(Ftinvvt);
    Zt.multVecTrans(1.0,rt,1.0,Ftinvvt);
    
    
    
    /*****************
    This calculates $$N_{t-1} = Z^T_tF_t^{-1}Z_t+L^T_tN_tL_t.$$
    *****************/
    GeneralMatrix NtLt(Nt,Lt);
    Nt.zeros();
    Nt.multAndAdd(Lt,"trans",NtLt);
    GeneralMatrix FtinvZt(Zt);
    Ftinv.multInvLeft(FtinvZt);
    Nt.multAndAdd(Zt,"trans",ConstGeneralMatrix(FtinvZt));
    
    
    
    /*****************
    This calculates $$\alpha_t = a_t + P_tr_{t-1}.$$
    *****************/
    alphat= (const Vector&)at;
    Pt.multVec(1.0,alphat,1.0,rt);
    
    
    /*****************
    This calculates $$V_t = P_t - P_tN_{t-1}P_t.$$
    *****************/
    Vt= (const GeneralMatrix&)Pt;
    GeneralMatrix NtPt(Nt,Pt);
    Vt.multAndAdd(Pt,NtPt,-1.0);
    
  }

  /*****************
  The non-diffuse smoother just performs a series of
  |KalmanTask::smootherNonDiffuseStep|.
  *****************/
  void 
  KalmanTask::smootherNonDiffuse(const FilterResults&fres,
    SmootherResults&sres)const
    {
    Vector rt(ssf.m);
    rt.zeros();
    GeneralMatrix Nt(ssf.m,ssf.m);
    Nt.zeros();
    for(int t= data.numCols();t>=1;t--)
      {
      Vector alphat(ssf.m);
      GeneralMatrix Vt(ssf.m,ssf.m);
      Vector etat(ssf.r);
      smootherNonDiffuseStep(t,fres,rt,Nt,alphat,Vt,etat);
      sres.set(t,alphat,etat,Vt);
      }
    }
  
  /*****************
  Here we cycle from $t=T,\ldots, 1$. Whenever $P_\infty$ is zero, we
  perform the non-diffuse step. Otherwise we permorn a common code to
  diffuse smoothing and then fork according to regularity of $F_\infty$.
  *****************/
  void 
  KalmanTask::smootherDiffuse(const DiffuseFilterResults&fres,
    SmootherResults&sres)const
    {
    Vector rt_0(ssf.m);
    Vector rt_1(ssf.m);
    GeneralMatrix Nt_0(ssf.m,ssf.m);
    GeneralMatrix Nt_1(ssf.m,ssf.m);
    GeneralMatrix Nt_2(ssf.m,ssf.m);
    rt_0.zeros();
    rt_1.zeros();
    Nt_0.zeros();
    Nt_1.zeros();
    Nt_2.zeros();
    
    for(int t= data.numCols();t>=1;t--)
      {
      Vector alphat(ssf.m);
      GeneralMatrix Vt(ssf.m,ssf.m);
      Vector etat(ssf.r);
      if(fres.isPinfZero(t))
        {
        smootherNonDiffuseStep(t,fres,rt_0,Nt_0,alphat,Vt,etat);
        }
      else
        {
        const Vector&vt= fres.getV(t);
        const GeneralMatrix&Lt_0= fres.getL0(t);
        const Vector&at= fres.getA(t);
        const GeneralMatrix&Ptstar= fres.getPstar(t);
        const GeneralMatrix&Ptinf= fres.getPinf(t);
        ConstGeneralMatrix Zt(((const TMatrix&)*ssf.Z)[t]);
        ConstGeneralMatrix Qt(((const TMatrix&)*ssf.Q)[t]);
        ConstGeneralMatrix Tt(((const TMatrix&)*ssf.T)[t]);
        ConstGeneralMatrix Rt(((const TMatrix&)*ssf.R)[t]);
        bool isTunit= ssf.T->isUnit(t);
        bool isQzero= ssf.Q->isZero(t);
        bool isRzero= ssf.R->isZero(t);
        
        /*****************
        Calculate $$\eta_t =Q_tR_tr_t^{(0)}.$$
        *****************/
        etat.zeros();
        if(!isQzero&&!isRzero)
          {
          Rt.multVecTrans(0.0,etat,1.0,rt_0);
          Vector etatsav((const Vector&)etat);
          Qt.multVec(0.0,etat,1.0,etatsav);
          }
        
        if(!fres.isFinfRegular(t))
          {
          /*****************
          We call here |smootherNonDiffuseStep| and calculate $r_{t-1}^{(1)}$,
          $N_{t-1}^{(1)}$, $N_{t-1}^{(2)}$ and correct for $\alpha_t$.
          *****************/
          smootherNonDiffuseStep(t,fres,rt_0,Nt_0,alphat,Vt,etat);
          
          /*****************
          This calculates $$r_{t-1}^{(1)} = T^T_tr_t^{(1)}.$$
          *****************/
          if(!isTunit)
            {
            Vector rt_1sav((const Vector&)rt_1);
            rt_1.zeros();
            Tt.multVecTrans(0.0,rt_1,1.0,rt_1sav);
            }
          
            /*****************
            This corrects $\alpha_t$ after|KalmanTask::smootherNonDiffuseStep|. 
            This adds $P_{\infty,t}r_{t-1}^{(1)}$ to $\alpha_t$.
          *****************/
          Ptinf.multVec(1.0,alphat,1.0,rt_1);
          
          /*****************
          This calculates $$N_{t-1}^{(1)} = T_t^TN_t^{(1)}L_t^{(0)}.$$
          *****************/
          if(!isTunit)
            {
            GeneralMatrix Nt_1Lt_0(Nt_1,Lt_0);
            Nt_1.zeros();
            Nt_1.multAndAdd(Tt,"trans",ConstGeneralMatrix(Nt_1Lt_0));
            }
          else
            Nt_1.mult(Nt_1,Lt_0);
          
            /*****************
            This calculates $$N_{t-1}^{(2)} = T_t^TN_t^{(2)}T_t.$$
          *****************/
          if(!isTunit)
            {
            GeneralMatrix Nt_2Tt(Nt_2,Tt);
            Nt_2.zeros();
            Nt_2.multAndAdd(Tt,"trans",ConstGeneralMatrix(Nt_2Tt));
            }
          
          }
        else
          {
          const GeneralMatrix&Lt_1= fres.getL1(t);
          const GeneralMatrix&Ft_2= fres.getF2(t);
          const PLUFact&Ftinfinv= fres.getFinfInverse(t);
          
          /*****************
          This calculates $$r_{t-1}^{(1)} = Z^T_tF_{\infty,t}^{-1}v_t^{(0)} +
          L_t^{(0)T}r_t^{(1)} + L_t^{(1)T}r_t^{(0)}.$$
          *****************/
          Vector rt_1sav((const Vector&)rt_1);
          Lt_0.multVecTrans(0.0,rt_1,1.0,rt_1sav);
          Lt_1.multVecTrans(1.0,rt_1,1.0,rt_0);
          Vector Ftinfinvvt(vt);
          Ftinfinv.multInvLeft(Ftinfinvvt);
          Zt.multVecTrans(1.0,rt_1,1.0,Ftinfinvvt);
          
          /*****************
          This calculates $$r_{t-1}^{(0)} = L_t^{(0)}r_t^{(0)}.$$
          *****************/
          Vector rt_0sav((const Vector&)rt_0);
          Lt_0.multVecTrans(0.0,rt_0,1.0,rt_0sav);
          
          /*****************
          This calculates 
          $$N_{t-1}^{(2)} = Z_t^TF_t^{(2)}Z_t + L_t^{(0)T}N_t^{(2)}L_t^{(0)}+
          L_t^{(0)T}N_t^{(1)}L_t^{(1)} + L_t^{(1)T}N_t^{(1)}L_t^{(0)}
          + L_t^{(1)T}N_t^{(0)}L_t^{(1)}.
          $$        
          *****************/
          GeneralMatrix Nt_2sav((const GeneralMatrix&)Nt_2);
          Nt_2.zeros();
          GeneralMatrix Ft_2Zt(Ft_2,Zt);
          Nt_2.multAndAdd(Zt,"trans",ConstGeneralMatrix(Ft_2Zt));
          GeneralMatrix Nt_2Lt_0(Nt_2sav,Lt_0);
          Nt_2.multAndAdd(Lt_0,"trans",Nt_2Lt_0);
          GeneralMatrix Nt_1Lt_1(Nt_1,Lt_1);
          Nt_2.multAndAdd(Lt_0,"trans",Nt_1Lt_1);
          GeneralMatrix Nt_1Lt_0(Nt_1,Lt_0);
          Nt_2.multAndAdd(Lt_1,"trans",Nt_1Lt_0);
          GeneralMatrix Nt_0Lt_1(Nt_0,Lt_1);
          Nt_2.multAndAdd(Lt_1,"trans",Nt_0Lt_1);
          
          /*****************
          This calculates 
          $$N_{t-1}^{(1)} = Z_t^TF_t^{(1)}Z_t + L_t^{(0)T}N_t^{(1)}L_t^{(0)}+
          L_t^{(1)T}N_t^{(0)}L_t^{(0)}.$$ |Nt_1Lt_0| was set in |@<calculate
          $N_{t-1}^{(2)}$ for diffuse smoother and regular $F_{\infty,t}$@>|.
          *****************/
          Nt_1.zeros();
          GeneralMatrix FtinfinvZt(Zt);
          Ftinfinv.multInvLeft(FtinfinvZt);
          Nt_1.multAndAdd(Zt,"trans",ConstGeneralMatrix(FtinfinvZt));
          Nt_1.multAndAdd(Lt_0,"trans",Nt_1Lt_0);
          GeneralMatrix Nt_0Lt_0(Nt_0,Lt_0);
          Nt_1.multAndAdd(Lt_1,"trans",Nt_0Lt_0);
          
          /*****************
          This calculates $$N_{t-1}^{(0)} = L_t^{(0)T}N_t^{(0)}L_t^{(0)}.$$
          |Nt_0Lt_0| was set in |@<calculate $N_{t-1}^{(1)}$ for diffuse
          smoother and regular $F_{\infty,t}$@>|.
          *****************/
          Nt_0.zeros();
          Nt_0.multAndAdd(Lt_0,"trans",Nt_0Lt_0);
          
          /*****************
          This calculates $$\alpha_t = a_t^{(0)} + P_{*,t}r_{t-1}^{(0)} + P_{\infty,t}r_{t-1}^{(1)}.$$
          for diffuse smoother and regular $F_{\infty,t}
          *****************/
          alphat= (const Vector&)at;
          Ptstar.multVec(1.0,alphat,1.0,rt_0);
          Ptinf.multVec(1.0,alphat,1.0,rt_1);
          }
        
          /*****************
          This calculates $$V_t = P_{*,t} - P_{*,t}N_{t-1}^{(0)}P_{*,t}
          - P_{\infty,t}N_{t-1}^{(1)}P_{*,t} -(P_{\infty,t}N_{t-1}^{(1)}P_{*,t})^T
          - P_{\infty,t}N_{t-1}^{(2)}P_{\infty,t}.$$
        *****************/
        Vt= (const GeneralMatrix&)Ptstar;
        GeneralMatrix Nt_0Ptstar(Nt_0,Ptstar);
        Vt.multAndAdd(Ptstar,Nt_0Ptstar,-1.0);
        GeneralMatrix Nt_2Ptinf(Nt_2,Ptinf);
        Vt.multAndAdd(Ptinf,Nt_2Ptinf,-1.0);
        GeneralMatrix Nt_1Ptstar(Nt_1,Ptstar);
        GeneralMatrix PtinfNt_1Ptstar(Ptinf,Nt_1Ptstar);
        Vt.add(-1.0,PtinfNt_1Ptstar);
        Vt.add(-1.0,PtinfNt_1Ptstar,"trans");
        
      }// end if/else
      sres.set(t,alphat,etat,Vt);
    }// end for
  }
  
  
  /*****************
  This evaluates a step loglikelihood as
  \log p(y_t\vert Y_{t-1})=-{1\over 2}\left[p\log(2\pi)+\log\vert F_t\vert+
  v_t^TF_t^{-1}v_t\right]$$
  This is a static method.
  *****************/
  
  double 
    KalmanTask::calcStepLogLik(const PLUFact&Finv,const Vector&v)
    {
    int p= Finv.numRows();
    Vector Finvv(v);
    Finv.multInvLeft(Finvv);
    double vFinvv= v.dot(Finvv);
    return-0.5*(p*log(2*M_PI)+Finv.getLogDeterminant()+vFinvv);
    }
  
  
  FilterUniResults::~FilterUniResults()
    {
    for(unsigned int i= 0;i<F.size();i++){
      if(L[i])
        delete L[i];
      if(a[i])
        delete a[i];
      if(P[i])
        delete P[i];
      }
    }
  
  ;
  
  void 
    FilterUniResults::set(int t,double FF,double vv,
    const GeneralMatrix&LL,const Vector&aa,
    const GeneralMatrix&PP,double ll)
    {
    TS_RAISE_IF(t<1||t> (int)L.size()+1,
      "Wrong time for FilterUniResults::set");
    
    int tm= t-1;
    if(L[tm])
      delete L[tm];
    if(a[tm])
      delete a[tm];
    if(P[tm])
      delete P[tm];
    
    if(t> maxt)
      maxt= t;
    
    F[tm]= FF;
    v[tm]= vv;
    L[tm]= new GeneralMatrix(LL);
    a[tm]= new Vector(aa);
    P[tm]= new GeneralMatrix(PP);
    loglik[tm]= ll;
    }
  
  ;
  
  double FilterUniResults::getF(int t)const
    {
    TS_RAISE_IF(t<1||t> maxt,
      "Wrong time for FilterUniResults::getF");
    return F[t-1];
    }
  
  ;
  
  double FilterUniResults::getV(int t)const
    {
    TS_RAISE_IF(t<1||t> maxt,
      "Wrong time for FilterUniResults::getV");
    return v[t-1];
    }
  
  ;
  
  const GeneralMatrix&FilterUniResults::getL(int t)const
    {
    TS_RAISE_IF(t<1||t> maxt,
      "Wrong time for FilterUniResults::getL");
    return*(L[t-1]);
    }
  
  ;
  
  const Vector&FilterUniResults::getA(int t)const
    {
    TS_RAISE_IF(t<1||t> maxt,
      "Wrong time for FilterUniResults::getA");
    return*(a[t-1]);
    }
  
  ;
  
  const GeneralMatrix&FilterUniResults::getP(int t)const
    {
    TS_RAISE_IF(t<1||t> maxt,
      "Wrong time for FilterUniResults::getP");
    return*(P[t-1]);
    }
  
    /*****************
    This adds all the log likelihoods for all periods. If some periods
    in the results have not been set, these are initialized to zeros and
    thus this method is pretty safe but may not be if the likelihood tends to
    be  far lower or higher than 0.
  *****************/
  double 
    FilterUniResults::getLogLikelihood()const
    {
    double res= 0.0;
    for(unsigned int i= 0;i<loglik.size();i++)
      res+= loglik[i];
    return res;
    }
  
  double 
    FilterUniResults::getLogLikelihood(int start)const
    {
    double res= 0.0;
    for(unsigned int i= start;i<loglik.size();i++)
      res+= loglik[i];
    return res;
    }
  
  double 
    FilterUniResults::getLogLikelihood(std::vector<double>* vloglik)const
    {
    double res= 0.0;
    for(unsigned int i= 0;i<loglik.size();i++)
      res+= loglik[i];
    *vloglik= loglik;
    return res;
    }
  
  double 
    FilterUniResults::getLogLikelihood(int start,std::vector<double>* vloglik)const
    {
    double res= 0.0;
    for(unsigned int i= start;i<loglik.size();i++)
      res+= loglik[i];
    *vloglik= loglik;
    return res;
    }
  
  
  DiffuseFilterUniResults::~DiffuseFilterUniResults()
    {
    for(unsigned int i= 0;i<L_1.size();i++){
      if(L_1[i])
        delete L_1[i];
      if(Pinf[i])
        delete Pinf[i];
      }
    }
  
  ;
  
  bool 
    DiffuseFilterUniResults::isFinfRegular(int t)const
    {
    TS_RAISE_IF(t<1||t> maxt,
      "Wrong time for DiffuseFilterUniResults::isFinfRegular");
    return Finf_reg[t-1];
    }
  
  ;
  
  bool 
    DiffuseFilterUniResults::isPinfZero(int t)const
    {
    TS_RAISE_IF(t<1||t> maxt,
      "Wrong time for DiffuseFilterUniResults::isPinfZero");
    return Pinf_zero[t-1];
    }
  
  ;
  
  double 
    DiffuseFilterUniResults::getFinf(int t)const
    {
    TS_RAISE_IF(t<1||t> maxt,
      "Wrong time for DiffuseFilterUniResults::getFinf");
    TS_RAISE_IF(!isFinfRegular(t),
      "Finf not regular in the period in DiffuseFilterUniResults::getFinf");
    return getF(t);
    }
  
  ;
  
  double 
    DiffuseFilterUniResults::getFstar(int t)const
    {
    TS_RAISE_IF(t<1||t> maxt,
      "Wrong time for DiffuseFilterUniResults::getFstar");
    TS_RAISE_IF(isFinfRegular(t),
      "Finf not zero in the period in DiffuseFilterUniResults::getFstar");
    return getF(t);
    }
  
  
  ;
  
  double 
    DiffuseFilterUniResults::getF2(int t)const
    {
    TS_RAISE_IF(t<1||t> maxt,
      "Wrong time for DiffuseFilterUniResults::getF2");
    TS_RAISE_IF(!isFinfRegular(t),
      "Finf not regular in the period in DiffuseFilterUniResults::getF2");
    return F_2[t-1];
    }
  
  ;
  
  const 
    GeneralMatrix&DiffuseFilterUniResults::getL1(int t)const
    {
    TS_RAISE_IF(t<1||t> maxt,
      "Wrong time for FilterUniResults::getL1");
    TS_RAISE_IF(!isFinfRegular(t),
      "Finf not regular in the period in DiffuseFilterUniResults::getL1");
    return*(L_1[t-1]);
    }
  
  ;
  
  const 
    GeneralMatrix&DiffuseFilterUniResults::getPinf(int t)const
    {
    TS_RAISE_IF(t<1||t> maxt,
      "Wrong time for FilterUniResults::getPinf");
    TS_RAISE_IF(isPinfZero(t),
      "Pinf is zero in the period in DiffuseFilterUniResults::getPinf");
    return*(Pinf[t-1]);
    }
  
  ;
  
  void 
    DiffuseFilterUniResults::set(int t,double FF,double FF_2,
    double vv,const GeneralMatrix&LL_0,
    const GeneralMatrix&LL_1,const Vector&aa,
    const GeneralMatrix&PPstar,const GeneralMatrix&PPinf,
    double ll)
    {
    FilterUniResults::set(t,FF,vv,LL_0,aa,PPstar,ll);
    
    int tm= t-1;
    if(L_1[tm])
      delete L_1[tm];
    if(Pinf[tm])
      delete Pinf[tm];
    
    L_1[tm]= new GeneralMatrix(LL_1);
    Pinf[tm]= new GeneralMatrix(PPinf);
    F_2[tm]= FF_2;
    Finf_reg[tm]= true;
    Pinf_zero[tm]= false;
    }
  
  ;
  
  void 
    DiffuseFilterUniResults::set(int t,double FFstar,double vv,
    const GeneralMatrix&LL_0,const Vector&aa,
    const GeneralMatrix&PPstar,const GeneralMatrix&PPinf,
    double ll)
    {
    FilterUniResults::set(t,FFstar,vv,LL_0,aa,PPstar,ll);
    
    int tm= t-1;
    if(Pinf[tm])
      delete Pinf[tm];
    Pinf[tm]= new GeneralMatrix(PPinf);
    
    Finf_reg[tm]= false;
    Pinf_zero[tm]= false;
    }
  
  ;
  
  int 
    DiffuseFilterUniResults::getDiffusePeriods()const
    {
    int d= maxt;
    while(d> 1&&isPinfZero(d))
      d--;
    return d;
    }
  
    /*****************
    @ This converts a multivariate |KalmanTask| to univariate
    |KalmanUniTask|. It unfolds time dimension so that at each time only
    one univariate observation comes. The measurment equation is
    transformed so that the measurment errors would not be correlated.
  *****************/
  KalmanUniTask::KalmanUniTask(const KalmanTask&kt)
    :me(kt.data,*(kt.ssf.Z),*(kt.ssf.H)),
    ssf(TMatrixCycle(*(me.Z),"rows"),TScalarCycle(*(me.H)),
    TMatrixPadUnit(*(kt.ssf.T),kt.data.numRows()),
    TMatrixPadZero(*(kt.ssf.R),kt.data.numRows()),
    TMatrixPadZero(*(kt.ssf.Q),kt.data.numRows())),
    data(me.y.base(),1,me.y.numRows()*me.y.numCols()),
    init(kt.init)
    {
    }
  
  
  double 
    KalmanUniTask::filter(int&per,int&d)const
    {
    if(!init.isDiffuse())
      {
      FilterUniResults fres(data.numCols());
      filterNonDiffuse(init.getA(),init.getPstar(),1,fres);
      d= 0;
      per= fres.getMaxT();
      if(fres.hasFinished())
        return fres.getLogLikelihood();
      }
    else
      {
      DiffuseFilterUniResults fres(data.numCols());
      filterDiffuse(init.getA(),init.getPstar(),init.getPinf(),
        1,fres);
      d= fres.getDiffusePeriods();
      per= fres.getMaxT();
      if(fres.hasFinished())
        return fres.getLogLikelihood();
      }
    return 0.0;
    }
  
  double 
    KalmanUniTask::filter(int&per,int&d, int start, std::vector<double>* vll)const
    {
    if(!init.isDiffuse())
      {
      FilterUniResults fres(data.numCols());
      filterNonDiffuse(init.getA(),init.getPstar(),1,fres);
      d= 0;
      per= fres.getMaxT();
      if(fres.hasFinished())
        return fres.getLogLikelihood(start,vll);
      }
    else
      {
      DiffuseFilterUniResults fres(data.numCols());
      filterDiffuse(init.getA(),init.getPstar(),init.getPinf(),
        1,fres);
      d= fres.getDiffusePeriods();
      per= fres.getMaxT();
      if(fres.hasFinished())
        return fres.getLogLikelihood(start,vll);
      }
    return 0.0;
    }
  
  double 
    KalmanUniTask::filter_and_smooth(SmootherResults&sres,
    int&per,int&d)const
    {
    if(!init.isDiffuse())
      {
      FilterUniResults fres(data.numCols());
      filterNonDiffuse(init.getA(),init.getPstar(),1,fres);
      d= 0;
      per= fres.getMaxT();
      if(fres.hasFinished())
        {
        smootherNonDiffuse(fres,sres);
        return fres.getLogLikelihood();
        }
      }
    else
      {
      DiffuseFilterUniResults fres(data.numCols());
      filterDiffuse(init.getA(),init.getPstar(),init.getPinf(),
        1,fres);
      d= fres.getDiffusePeriods();
      per= fres.getMaxT();
      if(fres.hasFinished())
        {
        smootherDiffuse(fres,sres);
        return fres.getLogLikelihood();
        }
      }
    return 0.0;
    }
  
    /*****************
    This filters univariate data starting at given $t$, $a_t$ and
    $P_t$. If at some period $F_t\leq 0$, than we end and the filter
    results are not finished.
  *****************/
  void 
    KalmanUniTask::filterNonDiffuse(const Vector&a,const GeneralMatrix&P,
    int first,FilterUniResults&fres)const
    {
    Vector at(a);
    GeneralMatrix Pt(P);
    for(int t= first;t<=data.numCols();t++){
      double yt= data.get(0,t-1);
      ConstGeneralMatrix Zt(((const TMatrix&)*ssf.Z)[t]);
      double Ht= ((const TScalar&)*ssf.H)[t];
      ConstGeneralMatrix Tt(((const TMatrix&)*ssf.T)[t]);
      ConstGeneralMatrix Qt(((const TMatrix&)*ssf.Q)[t]);
      ConstGeneralMatrix Rt(((const TMatrix&)*ssf.R)[t]);
      bool isTunit= ssf.T->isUnit(t);
      bool isQzero= ssf.Q->isZero(t);
      bool isRzero= ssf.R->isZero(t);
      
      
      double vt= at.dot(Zt.getData());
      vt= yt-vt;
      
      Vector Mt(ssf.m);
      Mt.zeros();
      Pt.multVec(0.0,Mt,1.0,Zt.getData());
      double Ft= Mt.dot(Zt.getData());
      Ft+= Ht;
      
      if(Ft<=0.0)
        return;
      
      Vector Kt(ssf.m);
      Kt.zeros();
      if(isTunit)
        Kt.add(1.0/Ft,Mt);
      else
        Tt.multVec(0.0,Kt,1.0/Ft,Mt);
      
      GeneralMatrix Lt(Tt);
      Lt.multAndAdd(ConstGeneralMatrix(Kt.base(),ssf.m,1),Zt,-1.0);
      
      double ll= calcStepLogLik(Ft,vt);
      fres.set(t,Ft,vt,Lt,at,Pt,ll);
      
      if(t<data.numCols())
        {
        
        if(!isTunit)
          {
          Vector atsave((const Vector&)at);
          Tt.multVec(0.0,at,1.0,atsave);
          }
        at.add(vt,Kt);
        
        
        GeneralMatrix PtLttrans(Pt,Lt,"trans");
        if(!isTunit)
          {
          Pt.zeros();
          Pt.multAndAdd(Tt,ConstGeneralMatrix(PtLttrans));
          }
        else
          {
          Pt= (const GeneralMatrix&)PtLttrans;
          }
        if(!isRzero&&!isQzero)
          {
          GeneralMatrix QtRttrans(Qt,Rt,"trans");
          Pt.multAndAdd(Rt,ConstGeneralMatrix(QtRttrans));
          }
        }
      }
    }
  
    /*****************
    This is a univariate version of |KalmanTask::filterDiffuse|. The
    decision whether $F_\infty$ is regular or zero is simpler here, yet
    the algorithm is not numerically stable. We recognize $d$ in the same
    way as in |KalmanTask::filterDiffuse|. It may still happen that small
    non-zero $F_t$ implies a wrong $P_{\infty,t+1}$, or zero $F_t$ which
    should be positive causes $d$ to be missed and numerical error is
    committed in $P_{\infty,t+1}$.
    
      So, as in |KalmanTask::filterDiffuse| we use |ndiff| and cancel it by
      one for periods for which $F_t$ is non-zero. If $P_\infty$ is not
      positive definite, we set it to zero.
  *****************/
  void 
    KalmanUniTask::filterDiffuse(const Vector&a,const GeneralMatrix&Pstar,
    const GeneralMatrix&Pinf,int first,
    DiffuseFilterUniResults&fres)const
    {
    Vector at(a);
    GeneralMatrix Ptstar(Pstar);
    GeneralMatrix Ptinf(Pinf);
    int ndiff= init.getNDiff();
    for(int t= first;t<=data.numCols();t++){
      
    /*****************
    This is the same code as |@<run non-diffuse filter 
    but it is semantically different, so we copy the code here.
      *****************/
      if(TSUtils::isZero(Ptinf))
        {
        filterNonDiffuse(at,Ptstar,t,fres);
        return;
        }
      
      double yt= data.get(0,t-1);
      ConstGeneralMatrix Zt(((const TMatrix&)*ssf.Z)[t]);
      double Ht= ((const TScalar&)*ssf.H)[t];
      ConstGeneralMatrix Tt(((const TMatrix&)*ssf.T)[t]);
      ConstGeneralMatrix Qt(((const TMatrix&)*ssf.Q)[t]);
      ConstGeneralMatrix Rt(((const TMatrix&)*ssf.R)[t]);
      bool isTunit= ssf.T->isUnit(t);
      bool isQzero= ssf.Q->isZero(t);
      bool isRzero= ssf.R->isZero(t);
      
      double vt= at.dot(Zt.getData());
      vt= yt-vt;
      
      Vector Mtstar(ssf.m);
      Mtstar.zeros();
      Ptstar.multVec(0.0,Mtstar,1.0,Zt.getData());
      
      double Ftstar= Mtstar.dot(Zt.getData());
      Ftstar+= Ht;
      
      Vector Mtinf(ssf.m);
      Mtinf.zeros();
      Ptinf.multVec(0.0,Mtinf,1.0,Zt.getData());
      
      double Ftinf= Mtinf.dot(Zt.getData());
      if(Ftinf<2*DBL_EPSILON)
        Ftinf= 0.0;
      
      if(Ftinf> 0.0)
        {
        ndiff--;
        
        double Ft_2= -Ftstar/Ftinf/Ftinf;
        
        Vector Kt_0(ssf.m);
        Kt_0.zeros();
        if(!isTunit)
          Tt.multVec(0.0,Kt_0,1.0/Ftinf,Mtinf);
        else
          Kt_0.add(1.0/Ftinf,Mtinf);
        
        Vector Kt_1tmp(ssf.m);
        Kt_1tmp.zeros();
        Kt_1tmp.add(Ft_2,Mtinf);
        Kt_1tmp.add(1.0/Ftinf,Mtstar);
        Vector Kt_1(ssf.m);
        if(!isTunit)
          {
          Kt_1.zeros();
          Tt.multVec(0.0,Kt_1,1.0,Kt_1tmp);
          }
        else
          {
          Kt_1= (const Vector&)Kt_1tmp;
          }
        
        GeneralMatrix Lt_0(Tt);
        Lt_0.multAndAdd(ConstGeneralMatrix(Kt_0.base(),ssf.m,1),Zt,-1.0);
        
        
        GeneralMatrix Lt_1(ConstGeneralMatrix(Kt_1.base(),ssf.m,1),Zt);
        Lt_1.mult(-1.0);
        
        
        double ll= -0.5*(log(2*M_PI)+log(Ftinf));
        fres.set(t,Ftinf,Ft_2,vt,Lt_0,Lt_1,at,Pstar,Pinf,ll);
        
        if(t<data.numCols())
          {
          /*****************
          This calculates $$a_{t+1} = T_ta_t+K_t^{(0)}v_t.$$
          *****************/
          if(!isTunit)
            {
            Vector atsave((const Vector&)at);
            Tt.multVec(0.0,at,1.0,atsave);
            }
          at.add(vt,Kt_0);
          
          GeneralMatrix tmp(Ptstar,Lt_0,"trans");
          tmp.multAndAdd(Ptinf,Lt_1,"trans");
          if(!isTunit)
            Ptstar.mult(Tt,ConstGeneralMatrix(tmp));
          else
            Ptstar= (const GeneralMatrix&)tmp;
          if(!isQzero&&!isRzero)
            {
            GeneralMatrix QtRttrans(Qt,Rt,"trans");
            Ptstar.multAndAdd(Rt,ConstGeneralMatrix(QtRttrans));
            }
          
          if(TSUtils::hasNegativeDiagonal(Ptstar)||!TSUtils::isSymDiagDominant(Ptstar))
            TSUtils::correctDefinitness(Ptstar);
          
          if(!isTunit)
            Ptinf.multLeft(Tt);
          Ptinf.multRightTrans(Lt_0);
          TSUtils::correctSymmetricity(Ptinf);
          
          
          
          if(TSUtils::hasNegativeDiagonal(Ptinf)||!TSUtils::isSymDiagDominant(Ptinf))
            TSUtils::correctDefinitness(Ptinf);
          
          if(ndiff==0)
            Ptinf.zeros();
          }
        
        }
      else
        {
        if(Ftstar==0.0)
          {
          return;
          }
        Vector Kt_0(ssf.m);
        Kt_0.zeros();
        if(!isTunit)
          Tt.multVec(0.0,Kt_0,1.0/Ftstar,Mtstar);
        else
          Kt_0.add(1.0/Ftstar,Mtstar);
        
        GeneralMatrix Lt_0(Tt);
        Lt_0.multAndAdd(ConstGeneralMatrix(Kt_0.base(),ssf.m,1),Zt,-1.0);
        
        double ll= calcStepLogLik(Ftstar,vt);
        fres.set(t,Ftstar,vt,Lt_0,at,Ptstar,Ptinf,ll);
        
        if(t<data.numCols())
          {
          if(!isTunit)
            {
            Vector atsave((const Vector&)at);
            Tt.multVec(0.0,at,1.0,atsave);
            }
          at.add(vt,Kt_0);
          
          if(!isTunit)
            {
            GeneralMatrix PtinfTttrans(Ptinf,Tt,"trans");
            Ptinf.mult(Tt,ConstGeneralMatrix(PtinfTttrans));
            }
          
          if(TSUtils::hasNegativeDiagonal(Ptinf)||!TSUtils::isSymDiagDominant(Ptinf))
            TSUtils::correctDefinitness(Ptinf);
          
          GeneralMatrix PtstarLt_0trans(Ptstar,Lt_0,"trans");
          if(!isTunit)
            Ptstar.mult(Tt,ConstGeneralMatrix(PtstarLt_0trans));
          else
            Ptstar= (const GeneralMatrix&)PtstarLt_0trans;
          if(!isQzero&&!isRzero)
            {
            GeneralMatrix QtRttrans(Qt,Rt,"trans");
            Ptstar.multAndAdd(Rt,ConstGeneralMatrix(QtRttrans));
            }
          
          
          if(TSUtils::hasNegativeDiagonal(Ptstar)||!TSUtils::isSymDiagDominant(Ptstar))
            TSUtils::correctDefinitness(Ptstar);
          
          }
        
        }
    }
  }
  
  
  void 
    KalmanUniTask::smootherNonDiffuseStep(int t,const FilterUniResults&fres,
    Vector&rt,GeneralMatrix&Nt,
    Vector&alphat,GeneralMatrix&Vt,
    Vector&etat)const
    {
    double Ft= fres.getF(t);
    double vt= fres.getV(t);
    const GeneralMatrix&Lt= fres.getL(t);
    const Vector&at= fres.getA(t);
    const GeneralMatrix&Pt= fres.getP(t);
    ConstGeneralMatrix Zt(((const TMatrix&)*ssf.Z)[t]);
    ConstGeneralMatrix Qt(((const TMatrix&)*ssf.Q)[t]);
    ConstGeneralMatrix Rt(((const TMatrix&)*ssf.R)[t]);
    bool isQzero= ssf.Q->isZero(t);
    bool isRzero= ssf.R->isZero(t);
    
    
    
    etat.zeros();
    if(!isQzero&&!isRzero)
      {
      Rt.multVecTrans(0.0,etat,1.0,rt);
      Vector etatsav((const Vector&)etat);
      Qt.multVec(0.0,etat,1.0,etatsav);
      }
    
    Vector rtsav((const Vector&)rt);
    Lt.multVecTrans(0.0,rt,1.0,rtsav);
    rt.add(vt/Ft,Zt.getData());
    
    GeneralMatrix NtLt(Nt,Lt);
    Nt.zeros();
    Nt.multAndAdd(Lt,"trans",NtLt);
    Nt.multAndAdd(Zt,"trans",Zt,1.0/Ft);
    
    alphat= (const Vector&)at;
    Pt.multVec(1.0,alphat,1.0,rt);
    
    Vt= (const GeneralMatrix&)Pt;
    GeneralMatrix NtPt(Nt,Pt);
    Vt.multAndAdd(Pt,NtPt,-1.0);
    }
  
  
  void
    KalmanUniTask::smootherNonDiffuse(const FilterUniResults&fres,
    SmootherResults&sres)const
    {
    Vector rt(ssf.m);
    rt.zeros();
    GeneralMatrix Nt(ssf.m,ssf.m);
    Nt.zeros();
    for(int t= data.numCols();t>=1;t--){
      Vector alphat(ssf.m);
      GeneralMatrix Vt(ssf.m,ssf.m);
      Vector etat(ssf.r);
      smootherNonDiffuseStep(t,fres,rt,Nt,alphat,Vt,etat);
      sres.set(t,alphat,etat,Vt);
      }
    }
  
  void 
    KalmanUniTask::smootherDiffuse(const DiffuseFilterUniResults&fres,
    SmootherResults&sres)const
    {
    
    Vector rt_0(ssf.m);
    Vector rt_1(ssf.m);
    GeneralMatrix Nt_0(ssf.m,ssf.m);
    GeneralMatrix Nt_1(ssf.m,ssf.m);
    GeneralMatrix Nt_2(ssf.m,ssf.m);
    rt_0.zeros();
    rt_1.zeros();
    Nt_0.zeros();
    Nt_1.zeros();
    Nt_2.zeros();
    
    
    for(int t= data.numCols();t>=1;t--)
      {
      Vector alphat(ssf.m);
      GeneralMatrix Vt(ssf.m,ssf.m);
      Vector etat(ssf.r);
      if(fres.isPinfZero(t))
        {
        smootherNonDiffuseStep(t,fres,rt_0,Nt_0,alphat,Vt,etat);
        }
      else
        {
        double vt= fres.getV(t);
        const GeneralMatrix&Lt_0= fres.getL0(t);
        const Vector&at= fres.getA(t);
        const GeneralMatrix&Ptstar= fres.getPstar(t);
        const GeneralMatrix&Ptinf= fres.getPinf(t);
        ConstGeneralMatrix Zt(((const TMatrix&)*ssf.Z)[t]);
        ConstGeneralMatrix Qt(((const TMatrix&)*ssf.Q)[t]);
        ConstGeneralMatrix Tt(((const TMatrix&)*ssf.T)[t]);
        ConstGeneralMatrix Rt(((const TMatrix&)*ssf.R)[t]);
        bool isTunit= ssf.T->isUnit(t);
        bool isQzero= ssf.Q->isZero(t);
        bool isRzero= ssf.R->isZero(t);
        
        etat.zeros();
        if(!isQzero&&!isRzero)
          {
          Rt.multVecTrans(0.0,etat,1.0,rt_0);
          Vector etatsav((const Vector&)etat);
          Qt.multVec(0.0,etat,1.0,etatsav);
          }
        
        if(!fres.isFinfRegular(t))
          {
          smootherNonDiffuseStep(t,fres,rt_0,Nt_0,alphat,Vt,etat);
          if(!isTunit)
            {
            Vector rt_1sav((const Vector&)rt_1);
            rt_1.zeros();
            Tt.multVecTrans(0.0,rt_1,1.0,rt_1sav);
            }
          
          Ptinf.multVec(1.0,alphat,1.0,rt_1);
          
          if(!isTunit)
            {
            GeneralMatrix Nt_1Lt_0(Nt_1,Lt_0);
            Nt_1.zeros();
            Nt_1.multAndAdd(Tt,"trans",ConstGeneralMatrix(Nt_1Lt_0));
            }
          else
            Nt_1.mult(Nt_1,Lt_0);
          
          
          if(!isTunit)
            {
            GeneralMatrix Nt_2Tt(Nt_2,Tt);
            Nt_2.zeros();
            Nt_2.multAndAdd(Tt,"trans",ConstGeneralMatrix(Nt_2Tt));
            }
          
          }
        else
          {
          const GeneralMatrix&Lt_1= fres.getL1(t);
          double Ft_2= fres.getF2(t);
          double Ftinf= fres.getFinf(t);
          
          
          Vector rt_1sav((const Vector&)rt_1);
          Lt_0.multVecTrans(0.0,rt_1,1.0,rt_1sav);
          Lt_1.multVecTrans(1.0,rt_1,1.0,rt_0);
          rt_1.add(vt/Ftinf,Zt.getData());
          
          
          Vector rt_0sav((const Vector&)rt_0);
          Lt_0.multVecTrans(0.0,rt_0,1.0,rt_0sav);
          
          GeneralMatrix Nt_2sav((const GeneralMatrix&)Nt_2);
          Nt_2.zeros();
          Nt_2.multAndAdd(Zt,"trans",Zt,Ft_2);
          GeneralMatrix Nt_2Lt_0(Nt_2sav,Lt_0);
          Nt_2.multAndAdd(Lt_0,"trans",Nt_2Lt_0);
          GeneralMatrix Nt_1Lt_1(Nt_1,Lt_1);
          Nt_2.multAndAdd(Lt_0,"trans",Nt_1Lt_1);
          GeneralMatrix Nt_1Lt_0(Nt_1,Lt_0);
          Nt_2.multAndAdd(Lt_1,"trans",Nt_1Lt_0);
          GeneralMatrix Nt_0Lt_1(Nt_0,Lt_1);
          Nt_2.multAndAdd(Lt_1,"trans",Nt_0Lt_1);
          
          
          
          Nt_1.zeros();
          Nt_1.multAndAdd(Zt,"trans",Zt,1.0/Ftinf);
          Nt_1.multAndAdd(Lt_0,"trans",Nt_1Lt_0);
          GeneralMatrix Nt_0Lt_0(Nt_0,Lt_0);
          Nt_1.multAndAdd(Lt_1,"trans",Nt_0Lt_0);
          
          Nt_0.zeros();
          Nt_0.multAndAdd(Lt_0,"trans",Nt_0Lt_0);
          
          alphat= (const Vector&)at;
          Ptstar.multVec(1.0,alphat,1.0,rt_0);
          Ptinf.multVec(1.0,alphat,1.0,rt_1);
          }
        
        Vt= (const GeneralMatrix&)Ptstar;
        GeneralMatrix Nt_0Ptstar(Nt_0,Ptstar);
        Vt.multAndAdd(Ptstar,Nt_0Ptstar,-1.0);
        GeneralMatrix Nt_2Ptinf(Nt_2,Ptinf);
        Vt.multAndAdd(Ptinf,Nt_2Ptinf,-1.0);
        GeneralMatrix Nt_1Ptstar(Nt_1,Ptstar);
        GeneralMatrix PtinfNt_1Ptstar(Ptinf,Nt_1Ptstar);
        Vt.add(-1.0,PtinfNt_1Ptstar);
        Vt.add(-1.0,PtinfNt_1Ptstar,"trans");
        
        
      }
      sres.set(t,alphat,etat,Vt);
    }
  }
  
  
  double 
    KalmanUniTask::calcStepLogLik(double F,double v)
    {
    return-0.5*(log(2*M_PI)+log(F)+v*v/F);
    }
  
  
  
