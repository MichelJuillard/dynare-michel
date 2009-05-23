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

#include "kalman.h"
#include "ts_exception.h"

#include <math.h> 
#include <float.h> 
#include <cmath> 



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


;

void FilterResults::set(int t,const PLUFact&FFinv,const Vector&vv,
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

;

const PLUFact&FilterResults::getFInverse(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getFInverse");
  return*(Finv[t-1]);
  }

;

const Vector&FilterResults::getV(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getV");
  return*(v[t-1]);
  }

;

const GeneralMatrix&FilterResults::getL(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getL");
  return*(L[t-1]);
  }

;

const Vector&FilterResults::getA(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getA");
  return*(a[t-1]);
  }

;

const GeneralMatrix&FilterResults::getP(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getP");
  return*(P[t-1]);
  }

;

double FilterResults::getLogLikelihood()const
  {
  double res= 0.0;
  for(unsigned int i= 0;i<loglik.size();i++)
    res+= loglik[i];
  return res;
  }

;

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

;

bool DiffuseFilterResults::isFinfRegular(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for DiffuseFilterResults::isFinfRegular");
  return Finf_reg[t-1];
  }

;

bool DiffuseFilterResults::isPinfZero(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for DiffuseFilterResults::isPinfZero");
  return Pinf_zero[t-1];
  }

;

const PLUFact&DiffuseFilterResults::getFinfInverse(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for DiffuseFilterResults::getFinfInverse");
  TS_RAISE_IF(!isFinfRegular(t),
    "Finf not regular in the period in DiffuseFilterResults::getFinfInverse");
  return getFInverse(t);
  }

;

const PLUFact&DiffuseFilterResults::getFstarInverse(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for DiffuseFilterResults::getFstarInverse");
  TS_RAISE_IF(isFinfRegular(t),
    "Finf not zero in the period in DiffuseFilterResults::getFstarInverse");
  return getFInverse(t);
  }

;

const GeneralMatrix&DiffuseFilterResults::getF2(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getF2");
  TS_RAISE_IF(!isFinfRegular(t),
    "Finf not regular in the period in DiffuseFilterResults::getF2");
  return*(F_2[t-1]);
  }

;

const GeneralMatrix&DiffuseFilterResults::getL1(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getL1");
  TS_RAISE_IF(!isFinfRegular(t),
    "Finf not regular in the period in DiffuseFilterResults::getL1");
  return*(L_1[t-1]);
  }

;

const GeneralMatrix&DiffuseFilterResults::getPinf(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterResults::getPinf");
  TS_RAISE_IF(isPinfZero(t),
    "Pinf is zero in the period in DiffuseFilterResults::getPinf");
  return*(L_1[t-1]);
  }

;

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

;

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

;

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


;

void SmootherResults::exportAlpha(GeneralMatrix&out)const
  {
  TS_RAISE_IF(mint> 1,
    "Results not finished in SmootherResults::exportAlpha");
  TS_RAISE_IF(out.numCols()!=(int)alpha.size(),
    "Wrong number of columns in SmootherResults::exportAlpha");
  TS_RAISE_IF(alpha[0]->length()!=out.numRows(),
    "Wrong number of rows in SmootherResults::exportAlpha");
  for(unsigned int tm= 0;tm<alpha.size();tm++){
    Vector outi(out,tm);
    outi= (const Vector&)(*alpha[tm]);
    }
  }

;

void SmootherResults::exportEta(GeneralMatrix&out)const
  {
  TS_RAISE_IF(mint> 1,
    "Results not finished in SmootherResults::exportEta");
  TS_RAISE_IF(out.numCols()!=(int)eta.size(),
    "Wrong number of columns in SmootherResults::exportEta");
  TS_RAISE_IF(eta[0]->length()!=out.numRows(),
    "Wrong number of rows in SmootherResults::exportEta");
  for(unsigned int tm= 0;tm<eta.size();tm++){
    Vector outi(out,tm);
    outi= (const Vector&)(*eta[tm]);
    }
  }

;

void SmootherResults::exportV(GeneralMatrix&out)const
  {
  TS_RAISE_IF(mint> 1,
    "Results not finished in SmootherResults::exportV");
  int m= V[0]->numRows();
  TS_RAISE_IF(out.numCols()!=(int)V.size()*m,
    "Wrong number of columns in SmootherResults::exportV");
  TS_RAISE_IF(m!=out.numRows(),
    "Wrong number of rows in SmootherResults::exportV");
  for(unsigned int tm= 0;tm<V.size();tm++){
    GeneralMatrix outi(out,0,tm*m,m,m);
    outi= (const GeneralMatrix&)(*V[tm]);
    }
  }


;

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

;

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

;

double KalmanTask::filter(int&per,int&d)const
  {
  if(!init.isDiffuse()){
    FilterResults fres(data.numCols());
    filterNonDiffuse(init.getA(),init.getPstar(),1,fres);
    d= 0;
    per= fres.getMaxT();
    if(fres.hasFinished())
      return fres.getLogLikelihood();
    }else{
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

;

double KalmanTask::filter_and_smooth(SmootherResults&sres,
                                     int&per,int&d)const
  {
  if(!init.isDiffuse()){
    FilterResults fres(data.numCols());
    filterNonDiffuse(init.getA(),init.getPstar(),1,fres);
    d= 0;
    per= fres.getMaxT();
    if(fres.hasFinished()){
      smootherNonDiffuse(fres,sres);
      return fres.getLogLikelihood();
      }
    }else{
    DiffuseFilterResults fres(data.numCols());
    filterDiffuse(init.getA(),init.getPstar(),init.getPinf(),
      1,fres);
    d= fres.getDiffusePeriods();
    per= fres.getMaxT();
    if(fres.hasFinished()){
      smootherDiffuse(fres,sres);
      return fres.getLogLikelihood();
      }
      }
    return 0.0;
  }

;

void KalmanTask::filterNonDiffuse(const Vector&a,const GeneralMatrix&P,
                                  int first,FilterResults&fres)const
  {
  Vector at(a);
  GeneralMatrix Pt(P);
  if(TSUtils::hasNegativeDiagonal(Pt)||!TSUtils::isSymDiagDominant(Pt))
    return;
  for(int t= first;t<=data.numCols();t++){
    ConstVector yt(data,t-1);
    ConstGeneralMatrix Zt(((const TMatrix&)*ssf.Z)[t]);
    ConstGeneralMatrix Ht(((const TMatrix&)*ssf.H)[t]);
    ConstGeneralMatrix Tt(((const TMatrix&)*ssf.T)[t]);
    ConstGeneralMatrix Qt(((const TMatrix&)*ssf.Q)[t]);
    ConstGeneralMatrix Rt(((const TMatrix&)*ssf.R)[t]);
    bool isTunit= ssf.T->isUnit(t);
    bool isQzero= ssf.Q->isZero(t);
    bool isRzero= ssf.R->isZero(t);
    
    
    Vector vt(yt);
    Zt.multsVec(vt,at);
    
    
    
    GeneralMatrix Mt(Pt,Zt,"trans");
    GeneralMatrix Ft(Ht);
    Ft.multAndAdd(Zt,ConstGeneralMatrix(Mt));
    
    
    PLUFact Ftinv(Ft);
    if(!Ftinv.isRegular())
      return;
    
    GeneralMatrix Kt(Tt,Mt);
    Ftinv.multInvRight(Kt);
    
    
    
    GeneralMatrix Lt(Tt);
    Lt.multAndAdd(ConstGeneralMatrix(Kt),Zt,-1.0);
    
    
    
    double ll= calcStepLogLik(Ftinv,vt);
    fres.set(t,Ftinv,vt,Lt,at,Pt,ll);
    
    
    if(t<data.numCols()){
      
      if(!isTunit){
        Vector atsave((const Vector&)at);
        Tt.multVec(0.0,at,1.0,atsave);
        }
      Kt.multVec(1.0,at,1.0,ConstVector(vt));
      
      
      
      GeneralMatrix PtLttrans(Pt,Lt,"trans");
      if(!isTunit){
        Pt.zeros();
        Pt.multAndAdd(Tt,ConstGeneralMatrix(PtLttrans));
        }else{
        Pt= (const GeneralMatrix&)PtLttrans;
          }
        if(!isRzero&&!isQzero){
          GeneralMatrix QtRttrans(Qt,Rt,"trans");
          Pt.multAndAdd(Rt,ConstGeneralMatrix(QtRttrans));
          }
        
        
      }
    }
  }

;

void KalmanTask::filterDiffuse(const Vector&a,const GeneralMatrix&Pstar,
                               const GeneralMatrix&Pinf,int first,
                               DiffuseFilterResults&fres)const
  {
  Vector at(a);
  GeneralMatrix Ptstar(Pstar);
  GeneralMatrix Ptinf(Pinf);
  int ndiff= init.getNDiff();
  for(int t= first;t<=data.numCols();t++){
    
    if(TSUtils::isZero(Ptinf)){
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
    
    
    Vector vt(yt);
    Zt.multsVec(vt,at);
    
    
    
    GeneralMatrix Mtstar(Ptstar,Zt,"trans");
    
    
    
    GeneralMatrix Ftstar(Ht);
    Ftstar.multAndAdd(Zt,ConstGeneralMatrix(Mtstar));
    
    
    
    GeneralMatrix Mtinf(Ptinf,Zt,"trans");
    
    
    
    GeneralMatrix Ftinf(Zt,ConstGeneralMatrix(Mtinf));
    
    
    PLUFact Ftinfinv(Ftinf);
    if(Ftinfinv.isRegular()&&Ftinfinv.getRcond()> 1.e-10){
      ndiff-= ssf.p;
      
      
      GeneralMatrix Ft_2(Ftstar);
      Ftinfinv.multInvRight(Ft_2);
      Ftinfinv.multInvLeft(Ft_2);
      Ft_2.mult(-1.0);
      
      
      
      GeneralMatrix Kt_0(Tt,Mtinf);
      Ftinfinv.multInvRight(Kt_0);
      
      
      
      GeneralMatrix Kt_1(Mtstar);
      Ftinfinv.multInvRight(Kt_1);
      Kt_1.multAndAdd(Mtinf,Ft_2);
      if(!isTunit)
        Kt_1.multLeft(Tt);
      
      
      
      GeneralMatrix Lt_0(Tt);
      Lt_0.multAndAdd(ConstGeneralMatrix(Kt_0),Zt,-1.0);
      
      
      
      GeneralMatrix Lt_1(Kt_1,Zt);
      Lt_1.mult(-1.0);
      
      
      
      double ll= -0.5*(ssf.p*log(2*M_PI)+Ftinfinv.getLogDeterminant());
      fres.set(t,Ftinfinv,Ft_2,vt,Lt_0,Lt_1,at,Pstar,Pinf,ll);
      
      
      if(t<data.numCols()){
        
        if(!isTunit){
          Vector atsave((const Vector&)at);
          Tt.multVec(0.0,at,1.0,atsave);
          }
        Kt_0.multVec(1.0,at,1.0,vt);
        
        
        
        GeneralMatrix tmp(Ptstar,Lt_0,"trans");
        tmp.multAndAdd(Ptinf,Lt_1,"trans");
        if(!isTunit)
          Ptstar.mult(Tt,ConstGeneralMatrix(tmp));
        else
          Ptstar= (const GeneralMatrix&)tmp;
        if(!isQzero&&!isRzero){
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
        
        
        
        if(ndiff<=0)
          Ptinf.zeros();
        }
      
      
      
      }else if(TSUtils::isZero(Ftinf)){
      
      PLUFact Ftstarinv(Ftstar);
      if(!Ftstarinv.isRegular()){
        return;
        }
      
      GeneralMatrix Kt_0(Tt,Mtstar);
      Ftstarinv.multInvRight(Kt_0);
      
      
      
      GeneralMatrix Lt_0(Tt);
      Lt_0.multAndAdd(ConstGeneralMatrix(Kt_0),Zt,-1.0);
      
      
      
      double ll= calcStepLogLik(Ftstarinv,vt);
      fres.set(t,Ftstarinv,vt,Lt_0,at,Ptstar,Ptinf,ll);
      
      
      if(t<data.numCols()){
        
        if(!isTunit){
          Vector atsave((const Vector&)at);
          Tt.multVec(0.0,at,1.0,atsave);
          }
        Kt_0.multVec(1.0,at,1.0,vt);
        
        
        
        if(!isTunit){
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
        if(!isQzero&&!isRzero){
          GeneralMatrix QtRttrans(Qt,Rt,"trans");
          Ptstar.multAndAdd(Rt,ConstGeneralMatrix(QtRttrans));
          }
        
        
        if(TSUtils::hasNegativeDiagonal(Ptstar)||!TSUtils::isSymDiagDominant(Ptstar))
          TSUtils::correctDefinitness(Ptstar);
        
        
        }
      
        }else{
        return;
          }
  }
}



void KalmanTask::smootherNonDiffuseStep(int t,const FilterResults&fres,
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
  
  
  etat.zeros();
  if(!isQzero&&!isRzero){
    Rt.multVecTrans(0.0,etat,1.0,rt);
    Vector etatsav((const Vector&)etat);
    Qt.multVec(0.0,etat,1.0,etatsav);
    }
  
  
  
  Vector rtsav((const Vector&)rt);
  Lt.multVecTrans(0.0,rt,1.0,rtsav);
  Vector Ftinvvt(vt);
  Ftinv.multInvLeft(Ftinvvt);
  Zt.multVecTrans(1.0,rt,1.0,Ftinvvt);
  
  
  
  GeneralMatrix NtLt(Nt,Lt);
  Nt.zeros();
  Nt.multAndAdd(Lt,"trans",NtLt);
  GeneralMatrix FtinvZt(Zt);
  Ftinv.multInvLeft(FtinvZt);
  Nt.multAndAdd(Zt,"trans",ConstGeneralMatrix(FtinvZt));
  
  
  
  alphat= (const Vector&)at;
  Pt.multVec(1.0,alphat,1.0,rt);
  
  
  Vt= (const GeneralMatrix&)Pt;
  GeneralMatrix NtPt(Nt,Pt);
  Vt.multAndAdd(Pt,NtPt,-1.0);
  
  
  
  }

;

void KalmanTask::smootherNonDiffuse(const FilterResults&fres,
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

;

void KalmanTask::smootherDiffuse(const DiffuseFilterResults&fres,
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
  
  
  for(int t= data.numCols();t>=1;t--){
    Vector alphat(ssf.m);
    GeneralMatrix Vt(ssf.m,ssf.m);
    Vector etat(ssf.r);
    if(fres.isPinfZero(t)){
      smootherNonDiffuseStep(t,fres,rt_0,Nt_0,alphat,Vt,etat);
      }else{
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
      
      
      etat.zeros();
      if(!isQzero&&!isRzero){
        Rt.multVecTrans(0.0,etat,1.0,rt_0);
        Vector etatsav((const Vector&)etat);
        Qt.multVec(0.0,etat,1.0,etatsav);
        }
      
      
      if(!fres.isFinfRegular(t)){
        
        smootherNonDiffuseStep(t,fres,rt_0,Nt_0,alphat,Vt,etat);
        
        if(!isTunit){
          Vector rt_1sav((const Vector&)rt_1);
          rt_1.zeros();
          Tt.multVecTrans(0.0,rt_1,1.0,rt_1sav);
          }
        
        
        
        Ptinf.multVec(1.0,alphat,1.0,rt_1);
        
        
        
        if(!isTunit){
          GeneralMatrix Nt_1Lt_0(Nt_1,Lt_0);
          Nt_1.zeros();
          Nt_1.multAndAdd(Tt,"trans",ConstGeneralMatrix(Nt_1Lt_0));
          }else
            Nt_1.mult(Nt_1,Lt_0);
          
          
          
          
          if(!isTunit){
            GeneralMatrix Nt_2Tt(Nt_2,Tt);
            Nt_2.zeros();
            Nt_2.multAndAdd(Tt,"trans",ConstGeneralMatrix(Nt_2Tt));
            }
          
          
          
          
          
        }else{
        
        const GeneralMatrix&Lt_1= fres.getL1(t);
        const GeneralMatrix&Ft_2= fres.getF2(t);
        const PLUFact&Ftinfinv= fres.getFinfInverse(t);
        
        
        Vector rt_1sav((const Vector&)rt_1);
        Lt_0.multVecTrans(0.0,rt_1,1.0,rt_1sav);
        Lt_1.multVecTrans(1.0,rt_1,1.0,rt_0);
        Vector Ftinfinvvt(vt);
        Ftinfinv.multInvLeft(Ftinfinvvt);
        Zt.multVecTrans(1.0,rt_1,1.0,Ftinfinvvt);
        
        
        
        Vector rt_0sav((const Vector&)rt_0);
        Lt_0.multVecTrans(0.0,rt_0,1.0,rt_0sav);
        
        
        
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
        
        
        
        Nt_1.zeros();
        GeneralMatrix FtinfinvZt(Zt);
        Ftinfinv.multInvLeft(FtinfinvZt);
        Nt_1.multAndAdd(Zt,"trans",ConstGeneralMatrix(FtinfinvZt));
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

;

double KalmanTask::calcStepLogLik(const PLUFact&Finv,const Vector&v)
  {
  int p= Finv.numRows();
  Vector Finvv(v);
  Finv.multInvLeft(Finvv);
  double vFinvv= v.dot(Finvv);
  return-0.5*(p*log(2*M_PI)+Finv.getLogDeterminant()+vFinvv);
  }

;

;


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

void FilterUniResults::set(int t,double FF,double vv,
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

;

double FilterUniResults::getLogLikelihood()const
  {
  double res= 0.0;
  for(unsigned int i= 0;i<loglik.size();i++)
    res+= loglik[i];
  return res;
  }


;

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

bool DiffuseFilterUniResults::isFinfRegular(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for DiffuseFilterUniResults::isFinfRegular");
  return Finf_reg[t-1];
  }

;

bool DiffuseFilterUniResults::isPinfZero(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for DiffuseFilterUniResults::isPinfZero");
  return Pinf_zero[t-1];
  }

;

double DiffuseFilterUniResults::getFinf(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for DiffuseFilterUniResults::getFinf");
  TS_RAISE_IF(!isFinfRegular(t),
    "Finf not regular in the period in DiffuseFilterUniResults::getFinf");
  return getF(t);
  }

;

double DiffuseFilterUniResults::getFstar(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for DiffuseFilterUniResults::getFstar");
  TS_RAISE_IF(isFinfRegular(t),
    "Finf not zero in the period in DiffuseFilterUniResults::getFstar");
  return getF(t);
  }


;

double DiffuseFilterUniResults::getF2(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for DiffuseFilterUniResults::getF2");
  TS_RAISE_IF(!isFinfRegular(t),
    "Finf not regular in the period in DiffuseFilterUniResults::getF2");
  return F_2[t-1];
  }

;

const GeneralMatrix&DiffuseFilterUniResults::getL1(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterUniResults::getL1");
  TS_RAISE_IF(!isFinfRegular(t),
    "Finf not regular in the period in DiffuseFilterUniResults::getL1");
  return*(L_1[t-1]);
  }

;

const GeneralMatrix&DiffuseFilterUniResults::getPinf(int t)const
  {
  TS_RAISE_IF(t<1||t> maxt,
    "Wrong time for FilterUniResults::getPinf");
  TS_RAISE_IF(isPinfZero(t),
    "Pinf is zero in the period in DiffuseFilterUniResults::getPinf");
  return*(Pinf[t-1]);
  }

;

void DiffuseFilterUniResults::set(int t,double FF,double FF_2,
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

void DiffuseFilterUniResults::set(int t,double FFstar,double vv,
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

int DiffuseFilterUniResults::getDiffusePeriods()const
  {
  int d= maxt;
  while(d> 1&&isPinfZero(d))
    d--;
  return d;
  }


;

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

;

double KalmanUniTask::filter(int&per,int&d)const
  {
  if(!init.isDiffuse()){
    FilterUniResults fres(data.numCols());
    filterNonDiffuse(init.getA(),init.getPstar(),1,fres);
    d= 0;
    per= fres.getMaxT();
    if(fres.hasFinished())
      return fres.getLogLikelihood();
    }else{
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

;

double KalmanUniTask::filter_and_smooth(SmootherResults&sres,
                                        int&per,int&d)const
  {
  if(!init.isDiffuse()){
    FilterUniResults fres(data.numCols());
    filterNonDiffuse(init.getA(),init.getPstar(),1,fres);
    d= 0;
    per= fres.getMaxT();
    if(fres.hasFinished()){
      smootherNonDiffuse(fres,sres);
      return fres.getLogLikelihood();
      }
    }else{
    DiffuseFilterUniResults fres(data.numCols());
    filterDiffuse(init.getA(),init.getPstar(),init.getPinf(),
      1,fres);
    d= fres.getDiffusePeriods();
    per= fres.getMaxT();
    if(fres.hasFinished()){
      smootherDiffuse(fres,sres);
      return fres.getLogLikelihood();
      }
      }
    return 0.0;
  }

;

void KalmanUniTask::filterNonDiffuse(const Vector&a,const GeneralMatrix&P,
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
    
    
    if(t<data.numCols()){
      
      if(!isTunit){
        Vector atsave((const Vector&)at);
        Tt.multVec(0.0,at,1.0,atsave);
        }
      at.add(vt,Kt);
      
      
      
      
      GeneralMatrix PtLttrans(Pt,Lt,"trans");
      if(!isTunit){
        Pt.zeros();
        Pt.multAndAdd(Tt,ConstGeneralMatrix(PtLttrans));
        }else{
        Pt= (const GeneralMatrix&)PtLttrans;
          }
        if(!isRzero&&!isQzero){
          GeneralMatrix QtRttrans(Qt,Rt,"trans");
          Pt.multAndAdd(Rt,ConstGeneralMatrix(QtRttrans));
          }
        
        
        
        
      }
    }
  }

;

void KalmanUniTask::filterDiffuse(const Vector&a,const GeneralMatrix&Pstar,
                                  const GeneralMatrix&Pinf,int first,
                                  DiffuseFilterUniResults&fres)const
  {
  Vector at(a);
  GeneralMatrix Ptstar(Pstar);
  GeneralMatrix Ptinf(Pinf);
  int ndiff= init.getNDiff();
  for(int t= first;t<=data.numCols();t++){
    
    if(TSUtils::isZero(Ptinf)){
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
    
    if(Ftinf> 0.0){
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
      if(!isTunit){
        Kt_1.zeros();
        Tt.multVec(0.0,Kt_1,1.0,Kt_1tmp);
        }else{
        Kt_1= (const Vector&)Kt_1tmp;
          }
        
        
        GeneralMatrix Lt_0(Tt);
        Lt_0.multAndAdd(ConstGeneralMatrix(Kt_0.base(),ssf.m,1),Zt,-1.0);
        
        
        GeneralMatrix Lt_1(ConstGeneralMatrix(Kt_1.base(),ssf.m,1),Zt);
        Lt_1.mult(-1.0);
        
        
        double ll= -0.5*(log(2*M_PI)+log(Ftinf));
        fres.set(t,Ftinf,Ft_2,vt,Lt_0,Lt_1,at,Pstar,Pinf,ll);
        
        if(t<data.numCols()){
          
          if(!isTunit){
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
          if(!isQzero&&!isRzero){
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
        
      }else{
      if(Ftstar==0.0){
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
      
      if(t<data.numCols()){
        
        if(!isTunit){
          Vector atsave((const Vector&)at);
          Tt.multVec(0.0,at,1.0,atsave);
          }
        at.add(vt,Kt_0);
        
        if(!isTunit){
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
        if(!isQzero&&!isRzero){
          GeneralMatrix QtRttrans(Qt,Rt,"trans");
          Ptstar.multAndAdd(Rt,ConstGeneralMatrix(QtRttrans));
          }
        
        
        if(TSUtils::hasNegativeDiagonal(Ptstar)||!TSUtils::isSymDiagDominant(Ptstar))
          TSUtils::correctDefinitness(Ptstar);
        
        
        }
      
        }
  }
}



void KalmanUniTask::smootherNonDiffuseStep(int t,const FilterUniResults&fres,
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
  if(!isQzero&&!isRzero){
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

;

void KalmanUniTask::smootherNonDiffuse(const FilterUniResults&fres,
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

;

void KalmanUniTask::smootherDiffuse(const DiffuseFilterUniResults&fres,
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
  
  
  for(int t= data.numCols();t>=1;t--){
    Vector alphat(ssf.m);
    GeneralMatrix Vt(ssf.m,ssf.m);
    Vector etat(ssf.r);
    if(fres.isPinfZero(t)){
      smootherNonDiffuseStep(t,fres,rt_0,Nt_0,alphat,Vt,etat);
      }else{
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
      if(!isQzero&&!isRzero){
        Rt.multVecTrans(0.0,etat,1.0,rt_0);
        Vector etatsav((const Vector&)etat);
        Qt.multVec(0.0,etat,1.0,etatsav);
        }
      
      
      
      
      if(!fres.isFinfRegular(t)){
        smootherNonDiffuseStep(t,fres,rt_0,Nt_0,alphat,Vt,etat);
        
        
        if(!isTunit){
          Vector rt_1sav((const Vector&)rt_1);
          rt_1.zeros();
          Tt.multVecTrans(0.0,rt_1,1.0,rt_1sav);
          }
        
        
        
        
        
        
        Ptinf.multVec(1.0,alphat,1.0,rt_1);
        
        
        
        
        
        
        if(!isTunit){
          GeneralMatrix Nt_1Lt_0(Nt_1,Lt_0);
          Nt_1.zeros();
          Nt_1.multAndAdd(Tt,"trans",ConstGeneralMatrix(Nt_1Lt_0));
          }else
            Nt_1.mult(Nt_1,Lt_0);
          
          
          
          
          
          
          
          if(!isTunit){
            GeneralMatrix Nt_2Tt(Nt_2,Tt);
            Nt_2.zeros();
            Nt_2.multAndAdd(Tt,"trans",ConstGeneralMatrix(Nt_2Tt));
            }
          
          
          
          
          
        }else{
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

;

double KalmanUniTask::calcStepLogLik(double F,double v)
  {
  return-0.5*(log(2*M_PI)+log(F)+v*v/F);
  }



