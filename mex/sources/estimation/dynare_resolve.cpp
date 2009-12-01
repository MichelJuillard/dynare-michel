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
/*****************************************************
% Based on Matlab Dynare
% function [int info] = dynare_resolve(ys,iv,ic,aux)
%
% function [A,B,ys,info] = dynare_resolve(iv,ic,aux)
% Computes the linear approximation and the matrices A and B of the
% transition equation and doing:
%	    check if ys is steady state
%	    dr
%     kalman_transition_matrix
%
% INPUTS
%    iv:             selected variables (observed and state variables)
%    ic:             state variables position in the transition matrix columns
%    aux:            indices for auxiliary equations
%
% MODIFIES
%    A:              matrix of predetermined variables effects in linear solution (ghx)
%    B:              matrix of shocks effects in linear solution (ghu)
%    ys:             steady state of original endogenous variables
%
% OUTPUTS
%    info=1:         the model doesn't determine the current variables '...' uniquely
%    info=2:         MJDGGES returns the following error code'
%    info=3:         Blanchard Kahn conditions are not satisfied: no stable '...' equilibrium
%    info=4:         Blanchard Kahn conditions are not satisfied:'...' indeterminacy
%    info=5:         Blanchard Kahn conditions are not satisfied:'...' indeterminacy due to rank failure
%    info=20:        can't find steady state info(2) contains sum of sqare residuals
%    info=30:        variance can't be computed
%
% SPECIAL REQUIREMENTS
%    none
********************************************************************/

#include "dsgeLikelihood.h"

int 
DsgeLikelihood::dynareResolveDR(vector<int>&iv,vector<int>&ic,GeneralMatrix& aux) // i.e. dynare_resolve()
  {
  //here comes a subset of [oo_.dr,info] = resol(oo_.steady_state,0);
  // check if ys is steady state and calcluate new one if not
  // testing for steadystate file: To Be Icluded at a later stage
  
#ifdef DEBUG
  mexPrintf("Calling SolveDRModel\n");
#endif 
  // [dr,info,M_,options_,oo_] = dr1(dr,check_flag,M_,options_,oo_);
  int DRinfo = SolveDRModel(endo_nbr,exo_nbr,nstatic, npred, nfwrd); //formerly known as dr1 in Matalb Dynare, i.e.
  if (DRinfo)
    {
    info = DRinfo;
    return DRinfo;
    }
  
  // End of resol: now rest of dynare_resolve:  
  
  // if nargin == 0
  //  if (iv.size()==0)
#ifdef DEBUG
  mexPrintf(" if iv==NULL\n");
#endif  
  if (&iv==NULL)
    {
#ifdef DEBUG
    mexPrintf(" construct iv\n");
#endif  
    //iv = (1:endo_nbr)';
    for (int i=1;i<=endo_nbr;++i) 
      iv.push_back(i);//= (1:endo_nbr)';
    }
  //  if (ic.size()==0)
#ifdef DEBUG
  mexPrintf(" if ic==NULL\n");
#endif  
  if (&ic==NULL)
    {
#ifdef DEBUG
    mexPrintf(" construct ic\n");
#endif  
    //ic = [ nstatic+(1:npred) endo_nbr+(1:size(oo_.dr.ghx,2)-npred) ]';
    //ic(npred+ghx.numCols());
    for(int i=0;i<npred;++i)
      {
      ic.push_back(i+nstatic+1);
      }
    for(int j=0;j<ghx.numCols()-npred;++j)
      ic.push_back(j+endo_nbr+1);
    }
  
#ifdef DEBUG
  mexPrintf(" if aux==NULL\n");
#endif  
  if (&aux==NULL)
    {
#ifdef DEBUG
    mexPrintf(" construct aux\n");
#endif  
    int i;
    aux =dr.getMatrixField(string("transition_auxiliary_variables"));
    //k = find(aux(:,2) > npred);
    //aux(:,2) = aux(:,2) + nstatic;
    vector<int>k(0);
    for (i=0; i< aux.numRows();++i)
      {
      if (aux.get(i,1)>npred) 
        k.push_back(i+1);
      aux.get(i, 1)+=nstatic;
      }
    
    //aux(k,2) = aux(k,2) + oo_.dr.nfwrd;
    for ( i=0; i< k.size();++i)
      aux.get(k[i]-1,1) +=nfwrd;
    }//end if
  
#ifdef DEBUG
  mexPrintf("[A,B] = kalman_transition_matrix\n");
#endif  
  
  // here is the content of [T R]=[A,B] = kalman_transition_matrix(oo_.dr,iv,ic,aux,M_.exo_nbr);
  
  int n_iv = iv.size();//length(iv);
  int n_ir1 = aux.numRows();// size(aux,1);
  int nr = n_iv + n_ir1;
  
  //    GeneralMatrix A=*(new GeneralMatrix (nr,nr));
  //    A.zeros();
  //    GeneralMatrix B=*(new GeneralMatrix(nr,exo_nbr));
  //    B.zeros();
  T.zeros();
  R.zeros();
  vector<int>i_n_iv(n_iv);
  for (int i=0;i<n_iv;++i) 
    i_n_iv[i]=i+1;//= (1:n_iv)';
  
#ifdef DEBUG
  mexPrintf("T=A assign by vec\n");
#endif  
  //A(i_n_iv,ic) = dr.ghx(iv,:);
  //A.AssignByVectors (i_n_iv,ic, ghx, iv, nullVec);//dr.ghx(iv,:);
  T.AssignByVectors (i_n_iv,ic, ghx, iv, nullVec);//dr.ghx(iv,:);
#ifdef DEBUG
  mexPrintf("Completed T=A assign by vec\n");
  T.print();
#endif    
  if (n_ir1 > 0)
    {
    //A(n_iv+1:end,:) = sparse(aux(:,1),aux(:,2),ones(n_ir1,1),n_ir1,nr);
    
#ifdef DEBUG
    mexPrintf("Create sparse\n");
#endif  
    GeneralMatrix sparse(n_ir1,nr);
    for (int i=0;i<n_ir1;++i)
      sparse.get((int)aux.get(i,0)-1,(int)aux.get(i,1)-1)=1;
    
      /* vector<int>span2end(A.numRows()-n_iv);
      for (int i=0;i<A.numRows()-n_iv;++i)
      span2end[i]=i+n_iv+1;
      A.place(sparse,n_iv,0); and T=A;
    */
#ifdef DEBUG
    mexPrintf("T.place (sparse,n_iv,0)\n");
#endif  
    T.place(sparse,n_iv,0);
    }
  
#ifdef DEBUG
  mexPrintf("R=B assign by vec R rows %d  cols %d  \n", R.numRows(), R.numCols());
#endif  
  //B(i_n_iv,:) = dr.ghu(iv,:); and R=B;
  R.AssignByVectors(i_n_iv, nullVec, ghu, iv, nullVec);
#ifdef DEBUG
  mexPrintf("Completed R=B assign by vec\n");
  R.print();
#endif
  
  //  ys = oo_.dr.ys;
  //    GeneralMatrix&ysmx = dr.getMatrixField(string("ys"));
  //    SteadyState=ysmx.getData();
  return DRinfo;
  }
