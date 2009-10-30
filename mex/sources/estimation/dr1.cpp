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
 
/***********************************************
% Based on Dynare Matlab
% function [dr,info,M_,options_,oo_] = SetDRModel(dr,task,M_,options_,oo_)
% formerly known as dr1 in Matalb Dynare
% Computes the reduced form solution of a rational expectation model (first or second order
% approximation of the stochastic model around the deterministic steady state). 
%
%  
********************************************************************/
#include "mexutils.h"
#include "DsgeLikelihood.h"

int 
DsgeLikelihood::SolveDRModel(const int endo_nbr, const int exo_nbr, const int nstatic, const int npred, int nfwrd)//dr1()
  {
  int  info = 0;
  int i;
 
  //    % expanding system for Optimal Linear Regulator
  if ((int)dynareParams.getDoubleField(string("ramsey_policy")))
    throw SYLV_MES_EXCEPTION("K-order-perturbation can not soleve for Ramsey policy");
  else     //    dr=set_state_space(dr,M_); - to be done prior to calling estaimte!!
    {    
    if (order ==1)
      {
      try
        {
        GeneralMatrix& ghx_u=SolveDRkOrderPert();//(dr,task,M_,options_, oo_ , ['.' mexext]);
        //SteadyState=ysteady;
        int sss= ghx_u.numCols();
        vector<int>span2nx(sss-exo_nbr); 
        for (i=0;i<sss-exo_nbr;++i)
          span2nx[i]=i+1;
        ghx= ( GeneralMatrix(ghx_u, nullVec,span2nx));//ghx_u(:,1:sss-M_.exo_nbr);  
        vector<int>spannx2end(exo_nbr); 
        for (i=0;i<exo_nbr;++i)
          spannx2end[i]=sss-exo_nbr+i+1;
        ghu= ( GeneralMatrix(ghx_u, nullVec,spannx2end)); //ghx_u(:,sss-M_.exo_nbr+1:end); 
        }
      catch(int e)
        {
        throw SYLV_MES_EXCEPTION("Problem with using k_order perturbation solver - Use Dynare solver instead");
        info = 4;
        penalty  = 1000; // info(2)
        return info;
        }//end
      }
    else //if options_.order >1
      {
      throw SYLV_MES_EXCEPTION("can not use order != 1 for estimation yet!");
      info = 4;
      penalty  = 1000;//info(2) = 1000;
      return info;
      };// end if
    
    if ((int)dynareParams.getDoubleField(string("loglinear")) == 1)
      {
      //k = find(dr.kstate(:,2) <= M_.maximum_endo_lag+1);
      vector<int>kk(0);
      int maximum_endo_lag=  (int)dynareParams.getDoubleField(string("maximum_endo_lag"));
      for(i=0;i<kstate.numRows();++i)
        if ( kstate.get(i,1)<=maximum_endo_lag+1)
          kk.push_back(i+1);
        
      //klag = dr.kstate(k,[1 2]);
      vector<int>kk2(2,1);
      kk2[1]=2;
      GeneralMatrix klag (kstate, kk,kk2);
      
      //k1 = dr.order_var;
      vector<int>k1klag(0);
      for (i=0; i< klag.numRows();++i)
        if ((int) klag.get(i, 0)>0)
          k1klag.push_back(order_var[(int) klag.get(i, 0)-1]);
      
      //dr.ghx = repmat(1./dr.ys(k1),1,size(dr.ghx,2)).*dr.ghx.* ...
      // ...repmat(dr.ys(k1(klag(:,1)))',size(dr.ghx,1),1);
      Vector invOrdSS(endo_nbr);//SteadyState.length());
      for (i=0;i<endo_nbr;++i)
        invOrdSS[i]=1/SteadyState[order_var[i-1]];
      GeneralMatrix mInvOrdSS(invOrdSS.base(),endo_nbr,1);
      GeneralMatrix&repInvSSx=mInvOrdSS.repmat(1,ghx.numCols());
      
      GeneralMatrix mSSt(SteadyState.base(),1,endo_nbr);
      
      GeneralMatrix mk1klagSSt(mSSt, k1klag,nullVec);
      GeneralMatrix&repk1klagSSt=mk1klagSSt.repmat(ghx.numRows(),1);
      
      ghx.multElements(repInvSSx);
      ghx.multElements(repk1klagSSt);
      
      //dr.ghu = repmat(1./dr.ys(k1),1,size(dr.ghu,2)).*dr.ghu;
      GeneralMatrix&repInvSSu=mInvOrdSS.repmat(1,ghu.numCols());
      ghu.multElements(repInvSSu);
      };//end if
    }
  }


/********************************************************************
* Solve the reduced DR k-order-model
*********************************************************************/
GeneralMatrix& 
DsgeLikelihood::SolveDRkOrderPert()  //kOrderPerturbation
  {
//  GeneralMatrix nullMat(0,0);
  model->getParams()=deepParams;
  model->getSteady()=SteadyState;
  
  try
    {
    
#ifdef DEBUG
    mexPrintf("k_order_perturbation: Calling walkStochSteady.\n");
#endif
    
    approx->walkStochSteady();
    
    //ConstTwoDMatrix ss(approx.getSS());
#ifdef DEBUG
    SteadyState=approx->getSS().getData();
    SteadyState.print();
#endif
    
    /* Write derivative outputs into memory map */
    map<string, ConstTwoDMatrix> mm;
    approx->getFoldDecisionRule().writeMMap(&mm);
    
#ifdef DEBUG
    approx->getFoldDecisionRule().print();
    mexPrintf("k_order_perturbation: Map print: \n");
    for (map<string, ConstTwoDMatrix>::const_iterator cit = mm.begin();
    cit != mm.end(); ++cit)
      {
      mexPrintf("k_order_perturbation: Map print: string: %s , g:\n", (*cit).first.c_str());
      (*cit).second.print();
      }
#endif
    // get latest ysteady
    SteadyState=model->getSteady();
#ifdef DEBUG
    SteadyState.print();
#endif
    
    // developement of the output.
#ifdef DEBUG
    mexPrintf("k_order_perturbation: Filling outputs.\n");
#endif
    int ii=1;
    GeneralMatrix* dgyu;
    /* Set the output pointer to the combined output matrix gyu = [gy gu]. */
    for (map<string, ConstTwoDMatrix>::const_iterator cit = mm.begin();
    ((cit != mm.end()) && (ii < 4)); ++cit)
      {
      if ((*cit).first!="g_0" && ii==2)
        {
        // TwoDMatrix dgyu((*cit).second.numRows(), (*cit).second.numCols(), mxGetPr(plhs[ii]));
#ifdef DEBUG
        dgyu=new GeneralMatrix ((*cit).second.numRows(), (*cit).second.numCols());
        *dgyu = (GeneralMatrix &)(*cit).second;
        mexPrintf("k_order_perturbation: cit %d print: \n", ii);
        (*cit).second.print();
        mexPrintf("k_order_perturbation: dguy %d print: \n", ii);
        dgyu->print(); //!! This print Crashes???
#endif
        return (GeneralMatrix &)(*cit).second;
        }
      ++ii;
      }
    return *dgyu;
    }
  catch (const KordException &e)
    {
    printf("Caugth Kord exception: ");
    e.print();
    mexPrintf("Caugth Kord exception: %s", e.get_message());
    }
  catch (const TLException &e)
    {
    printf("Caugth TL exception: ");
    e.print();
    }
  catch (SylvException &e)
    {
    printf("Caught Sylv exception: ");
    e.printMessage();
    }
  catch (const DynareException &e)
    {
    printf("Caught KordpDynare exception: %s\n", e.message());
    mexPrintf("Caugth Dynare exception: %s", e.message());
    }
  catch (const ogu::Exception &e)
    {
    printf("Caught ogu::Exception: ");
    e.print();
    mexPrintf("Caugth general exception: %s", e.message());
    }  //catch
  }; // end of mexFunction()
  
