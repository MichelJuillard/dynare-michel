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
  int  infoDR = 0;
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
        if (ghx_u.isZero()) 
          {
          mexPrintf("********  ghx_u is Zero ! *******\n");
          //throw(1);
          }
/**********/
        int sss= ghx_u.numCols();
#ifdef DEBUG
		mexPrintf("*********GHX_U colos %d then Allocate GHX and GHU *********\n",  sss);
		ghx_u.print();
#endif
        vector<int>span2nx(sss-exo_nbr); 
        for (i=0;i<sss-exo_nbr;++i)
          span2nx[i]=i+1;
        //ghx= ( (ghx_u, nullVec,span2nx));//ghx_u(:,1:sss-M_.exo_nbr);  
        GeneralMatrix gh_x (ghx_u, nullVec,span2nx);//ghx_u(:,1:sss-M_.exo_nbr);  
#ifdef DEBUG
		mexPrintf("*********GH_X*********\n");
		gh_x.print();
#endif 
         ghx= gh_x;
#ifdef DEBUG
		mexPrintf("*********GHX*********\n");
		ghx.print();
#endif 
        vector<int>spannx2end(exo_nbr); 
        for (i=0;i<exo_nbr;++i)
          spannx2end[i]=sss-exo_nbr+i+1;
        ghu= ( GeneralMatrix(ghx_u, nullVec,spannx2end)); //ghx_u(:,sss-M_.exo_nbr+1:end); 
/**********
  //Test only:
  ghu=dr.getMatrixField(string("ghu"));
  ghx=dr.getMatrixField(string("ghx"));
********/
#ifdef DEBUG
	mexPrintf("*********GHU*********\n");
  ghu.print();
//  ghx.print();
#endif 

  // end test
        delete &ghx_u;        
        }
      catch(int e)
        {
        throw SYLV_MES_EXCEPTION("Problem with using k_order perturbation solver ");
        info = 4;
        penalty  = 1000; // info(2)
        infoDR=info;
        return infoDR;
        }//end
      }
    else //if options_.order >1
      {
      throw SYLV_MES_EXCEPTION("can not use order != 1 for estimation yet!");
      info = 4;
      penalty  = 1000;//info(2) = 1000;
      infoDR=info;
      return infoDR;
      };// end if
    
    if ((int)dynareParams.getDoubleField(string("loglinear")) == 1)
      {
#ifdef DEBUG
        mexPrintf("setting SolveDRModel loglinear results\n");
#endif 
      //k = find(dr.kstate(:,2) <= M_.maximum_endo_lag+1);
      vector<int>kk(0);
      int maximum_endo_lag=  (int)dynareParams.getDoubleField(string("maximum_endo_lag"));
      for(i=0;i<kstate.numRows();++i)
        if ( kstate.get(i,1)<=maximum_endo_lag+1)
          kk.push_back(i+1);
        
      //klag = dr.kstate(k,[1 2]);
      vector<int>kk2(2);
      kk2[0]=1;
      kk2[1]=2;
#ifdef DEBUG
      mexPrintf("setting klag for loglinear results\n");
#endif
      GeneralMatrix klag (kstate, kk,kk2);
      
      //k1 = dr.order_var;
      vector<int>k1klag(0);
#ifdef DEBUG
        mexPrintf("setting k1klag for loglinear results\n");
#endif
      for (i=0; i< klag.numRows();++i)
        if ((int) klag.get(i, 0)>0)
          k1klag.push_back(order_var[(int) klag.get(i, 0)-1]);
      
      //dr.ghx = repmat(1./dr.ys(k1),1,size(dr.ghx,2)).*dr.ghx.* ...
      // ...repmat(dr.ys(k1(klag(:,1)))',size(dr.ghx,1),1);
      Vector invOrdSS(endo_nbr);//SteadyState.length());
      for (i=0;i<endo_nbr;++i)
        invOrdSS[i]=1/SteadyState[order_var[i]-1];

#ifdef DEBUG
        mexPrintf("setting mInvOrdSS for loglinear results\n");
#endif
      GeneralMatrix mInvOrdSS(invOrdSS.base(),endo_nbr,1);
#ifdef DEBUG
        mInvOrdSS.print();
        mexPrintf("SolveDRModel  Call repmat 1 for loglinear ghx results\n");
#endif       
      GeneralMatrix&repInvSSx=mInvOrdSS.repmat(1,ghx.numCols());

      Vector k1klagSS(k1klag.size());
      for (i=0;i<k1klag.size();++i)
        k1klagSS[i]=SteadyState[k1klag[i]-1];

     // GeneralMatrix mSSt(SteadyState.base(),1,endo_nbr);
     // GeneralMatrix mk1klagSSt(mSSt, k1klag,nullVec);
      
      GeneralMatrix mk1klagSSt(k1klagSS.base(), 1,k1klag.size());
#ifdef DEBUG
        mk1klagSSt.print();
        repInvSSx.print();
        mexPrintf("SolveDRModel  Call repmat 2 for loglinear ghx results\n");
#endif       
      GeneralMatrix&repk1klagSSt=mk1klagSSt.repmat(ghx.numRows(),1);
#ifdef DEBUG
        mexPrintf("Final setting SolveDRModel loglinear ghx results\n");
#endif       
      ghx.multElements(repInvSSx);
      ghx.multElements(repk1klagSSt);
      
      //dr.ghu = repmat(1./dr.ys(k1),1,size(dr.ghu,2)).*dr.ghu;
#ifdef DEBUG
        mexPrintf("SolveDRModel  Call repmat 1 for loglinear ghu results\n");
#endif       
      GeneralMatrix&repInvSSu=mInvOrdSS.repmat(1,ghu.numCols());
#ifdef DEBUG
        mexPrintf("Final setting SolveDRModel loglinear ghu results\n");
#endif       
      ghu.multElements(repInvSSu);
      delete &repk1klagSSt;
      delete &repInvSSu;
      delete &repInvSSx;
#ifdef DEBUG
	mexPrintf("Final loglinear ghu and ghx results\n");
  ghu.print();
  ghx.print();
#endif 
      };//end if
    }
  return infoDR;
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
  GeneralMatrix* dgyu=new GeneralMatrix (ghu.numRows(), ghx.numCols()+ghu.numCols());
  dgyu->zeros();
  
  try
    {
    
#ifdef DEBUG
    mexPrintf(" DeepParams:.\n");
    deepParams.print();
    mexPrintf(" Calling walkStochSteady with Params:.\n");
    model->getParams().print();
#endif
    
    approx->walkStochSteady();
#ifdef DEBUG
    mexPrintf("End of walkStochSteady - write map.\n");
#endif
    // Write derivative outputs into memory map 
    map<string, ConstTwoDMatrix> mm;
    approx->getFoldDecisionRule().writeMMap(&mm);
#ifdef DEBUG
    mexPrintf("k_order_perturbation: Map print: \n");
    approx->getFoldDecisionRule().print();
#endif    
#ifdef DEBUG
    mexPrintf("k_order_perturbation: Map print: \n");
    for (map<string, ConstTwoDMatrix>::const_iterator cit = mm.begin();
    cit != mm.end(); ++cit)
      {
      mexPrintf("k_order_perturbation: Map print: string: %s , g:\n", (*cit).first.c_str());
      (*cit).second.print();
      }
#endif
  // get latest ysteady
//    SteadyState=model->getSteady();
#ifdef DEBUG
//    mexPrintf("Steady State\n");
//    SteadyState.print();
#endif
    
    
    // developement of the output.
#ifdef DEBUG
    mexPrintf("k_order_perturbation: Filling outputs.\n");
#endif
    int ii=1;
    // Set the output pointer to the combined output matrix gyu = [gy gu]. 
    for (map<string, ConstTwoDMatrix>::const_iterator cit = mm.begin();
    ((cit != mm.end()) && (ii < 4)); ++cit)
      {
      if ((*cit).first!="g_0" && ii==2)
        {
        dgyu->getData() = (*cit).second.getData();
#ifdef DEBUG
        mexPrintf("k_order_perturbation: cit %d numRows %d numCols %d print: \n", ii, (*cit).second.numRows(), (*cit).second.numCols());
        (*cit).second.print();
        mexPrintf("k_order_perturbation: dguy output %d print: \n", ii);
        dgyu->print(); //!! This print Crashes???
#endif
        return *dgyu;
        }
      ++ii;
      }
    return *dgyu;
    }
  catch (const KordException &e)
    {
    printf("Caugth Kord exception in SolveDRkOrderPert: ");
    e.print();
    mexPrintf("Caugth Kord exception: %s", e.get_message());
    }
  catch (const TLException &e)
    {
    printf("Caugth TL exception in SolveDRkOrderPert: ");
    mexPrintf("Caugth TL exception in SolveDRkOrderPert: ");
    e.print();
    }
  catch (SylvException &e)
    {
    printf("Caught Sylv exception in SolveDRkOrderPert: ");
    mexPrintf("Caught Sylv exception in SolveDRkOrderPert: ");
    e.printMessage();
    }
  catch (const DynareException &e)
    {
    printf("Caught KordpDynare exception in SolveDRkOrderPert: %s\n", e.message());
    mexPrintf("Caugth Dynare exception in SolveDRkOrderPert: %s", e.message());
    }
  catch (const ogu::Exception &e)
    {
    printf("Caught ogu::Exception in SolveDRkOrderPert: ");
    e.print();
    mexPrintf("Caugth general exception inSolveDRkOrderPert: %s", e.message());
    }  //catch
  }; // end of mexFunction()
  
