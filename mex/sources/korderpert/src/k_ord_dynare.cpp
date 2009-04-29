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

// GP, based on work by O.Kamenik

#include <vector>
#include "first_order.h"
#include "k_ord_dynare.h"

#include "mex.h" 

#include "memory_file.h"


//#include "k_order_perturbation.h"
#ifndef DYNVERSION
#define DYNVERSION "unknown"
#endif

/**************************************************************************************/
/*       Dynare DynamicModel class                                                                 */
/**************************************************************************************/
class KordpJacobian;

KordpDynare::KordpDynare(const char** endo,  int num_endo,
			   const char** exo, int nexog, int npar, //const char** par,
   			   Vector* ysteady, TwoDMatrix* vcov, Vector* inParams, int nstat,
			   int npred, int nforw, int nboth, const int jcols, const int nsteps, int norder, //const char* modName,
			   Journal& jr, DynamicModelDLL& dynamicDLL, double sstol, const vector<int>* var_order, 
			   const TwoDMatrix * llincidence, double criterium )
	: nStat(nstat), nBoth(nboth), nPred(npred), nForw(nforw), nExog(nexog), nPar(npar),
	nYs(npred + nboth), nYss(nboth + nforw),nY(num_endo), nJcols(jcols), nSteps(nsteps), nOrder(norder), 
	journal(jr),  dynamicDLL(dynamicDLL), ySteady(ysteady), vCov(vcov), params (inParams), 
	md(1), dnl(NULL), denl(NULL), dsnl(NULL), ss_tol(sstol), varOrder( var_order), 
	ll_Incidence(llincidence), qz_criterium(criterium) 
{
#ifdef DEBUG		
   mexPrintf("k_ord_dynare Dynare constructor: ny=%d, order=%d, nPar=%d .\n", nY,nOrder,nPar);
	for (int i = 0; i < nY; i++) {
        mexPrintf("k_ord_dynare calling DynareNameList names[%d]= %s.\n", i, endo[i] );}
	for (int i = 0; i < nPar; i++) {
        mexPrintf("k_ord_perturbation: Params[%d]= %g.\n", i, (*params)[i]);  }
	for (int i = 0; i < nY; i++) {
        mexPrintf("k_ord_perturbation: ysteady[%d]= %g.\n", i, (*ySteady)[i]);  }
    mexPrintf("k_ord_dynare: dynare constructor, trying namelists.\n");
#endif		
	try{
		dnl = new DynareNameList(*this, endo);
		denl = new DynareExogNameList(*this, exo);
#ifdef DEBUG		
    mexPrintf("k_ord_dynare: dynare constructor, trying StateNamelist.\n");
#endif		
		dsnl = new DynareStateNameList(*this, *dnl, *denl);

		JacobianIndices= ReorderDynareJacobianIndices( varOrder);

	//	Initialise ModelDerivativeContainer(*this, this->md, nOrder);
    	for (int iord = 1; iord <= nOrder; iord++) {
    		FSSparseTensor* t = new FSSparseTensor(iord, nY+nYs+nYss+nExog,nY);
        	md.insert(t);
        }

	}
	catch (...){
        mexPrintf("k_ord_dynare: dynare constructor, error in StateNamelist construction.\n");
        throw DynareException(__FILE__, __LINE__, string("Could not construct Name Lists. \n"));
    }
}

KordpDynare::KordpDynare(const KordpDynare& dynare)
	: nStat(dynare.nStat), nBoth(dynare.nBoth),	nPred(dynare.nPred), 
	nForw(dynare.nForw), nExog(dynare.nExog),  nPar(dynare.nPar),
	nYs(dynare.nYs), nYss(dynare.nYss),nY(dynare.nY), nJcols(dynare.nJcols), 
	nSteps(dynare.nSteps), nOrder(dynare.nOrder), journal(dynare.journal),
	dynamicDLL(dynare.dynamicDLL), //modName(dynare.modName),
	ySteady(NULL), params(NULL), vCov(NULL), md(dynare.md), 
	dnl(NULL), denl(NULL), dsnl(NULL), ss_tol(dynare.ss_tol), 
	varOrder(dynare.varOrder), ll_Incidence(dynare.ll_Incidence),
	JacobianIndices(dynare.JacobianIndices), qz_criterium(dynare.qz_criterium)
{
	ySteady = new Vector(*(dynare.ySteady));
	params = new Vector(*(dynare.params));
	vCov = new TwoDMatrix(*(dynare.vCov));
	dnl = new DynareNameList(dynare);//(*this);
	denl = new DynareExogNameList(dynare);//(*this);
	dsnl = new DynareStateNameList(*this, *dnl, *denl);
}

KordpDynare::~KordpDynare()
{
	if (ySteady)
		delete ySteady;
	if (params)
		delete params;
	if (vCov)
		delete vCov;
	if (dnl)
		delete dnl;
	if (dsnl)
		delete dsnl;
	if (denl)
		delete denl;
}


/** This clears the container of model derivatives and initializes it
 * inserting empty sparse tensors up to the given order. */
ModelDerivativeContainer::ModelDerivativeContainer(const KordpDynare& model, 
		TensorContainer<FSSparseTensor>& mod_ders, int order): md(mod_ders)
{
	md.clear();
	for (int iord = 1; iord <= order; iord++) {
		FSSparseTensor* t = new FSSparseTensor(iord, model.ny()+model.nys()+model.nyss()+model.nexog(), model.ny());
		md.insert(t);
	}
}


void KordpDynare::solveDeterministicSteady(Vector& steady)
{
	JournalRecordPair pa(journal);
	pa << "Non-linear solver for deterministic steady state By-passed " << endrec;
	/*************************;  GP Dec 08 by-pass
	KordpVectorFunction dvf(*this);
	KordpJacobian dj(*this);
	ogu::NLSolver nls(dvf, dj, 500, ss_tol, journal);
	int iter;
	if (! nls.solve(*ySteady, iter))
		throw DynareException(__FILE__, __LINE__,
							  "Could not obtain convergence in non-linear solver");
	***************************/
}

// evaluate system at given y_t=y_{t+1}=y_{t-1}, and given shocks x_t
void KordpDynare::evaluateSystem(Vector& out, const Vector& yy, const Vector& xx)
{
	dynamicDLL.eval( yy,  xx, //int nb_row_x, 
				 params, //int it_, 
                 out, NULL, NULL);

}

// evaluate system at given y^*_{t-1}, y_t, y^{**}_{t+1} and at
// exogenous x_t, all three vectors yym, yy, and yyp have the
// respective lengths of y^*_{t-1}, y_t, y^{**}_{t+1}

void KordpDynare::evaluateSystem(Vector& out, const Vector& yym, const Vector& yy,
							const Vector& yyp, const Vector& xx)
{
#ifdef DEBUG
	mexPrintf("k_order_dynaare.cpp: Call eval in EvaluateSystem\n");
#endif
	dynamicDLL.eval( yy,  xx, //int nb_row_x, 
				 params, //int it_, 
                 out, NULL, NULL);
}
/************************************************
* this is main derivative calculation functin that indirectly calls dynamic.dll
* which performs actual calculation and reorders
***************************************************/ 
//void KordpDynare::calcDerivatives(const Vector& yy, const TwoDMatrix& xx)
void KordpDynare::calcDerivatives(const Vector& yy, const Vector& xx)
{
   // Hessian TwoDMatrix *g2;
	TwoDMatrix *g2=NULL;
    //Jacobian 
	TwoDMatrix * g1=new TwoDMatrix(nY,nJcols); // generate g1 for jacobian  
	g1->zeros();

	if ((nJcols != g1->ncols()) && ( nY != g1->nrows())) {
		mexPrintf("k_ord_dynare.cpp: Error in calcDerivatives: Created wrong jacobian");
		return;
	}

    // Hessian TwoDMatrix *g2;
	if (nOrder>1){
	//TwoDMatrix * 
		g2=new TwoDMatrix(nY,nJcols*nJcols); // generate g2 for Hessian  
		g2->zeros();
    }
	Vector& out= *(new Vector(nY));
	out.zeros();
	const Vector * llxYYp; // getting around the constantness
	if ((nJcols - nExog) > yy.length()){
		llxYYp=  (LLxSteady( yy));
	} else {
		llxYYp= &yy;
	}
	const Vector & llxYY=*(llxYYp);

#ifdef DEBUG
	mexPrintf("k_order_dynaare.cpp: Call eval in calcDerivatives\n");
#endif
	try {
			dynamicDLL.eval( llxYY,  xx, //int nb_row_x, 
					params, //int it_, 
					out, g1, g2);
//		}
	}
	catch (...){
			mexPrintf("k_ord_dynare.cpp: Error in dynamicDLL.eval in calcDerivatives");
			return;
	}
	if ((nJcols!=g1->ncols()) && ( nY != g1->nrows())) {
			mexPrintf("k_ord_dynare.cpp: Error in calcDerivatives: dynamicDLL.eval returned wrong jacobian");
			return;
	}
	//	ReorderCols(g1, JacobianIndices); and populate container
	populateDerivativesContainer(g1,1,JacobianIndices);
	if (nOrder>1){
	  //		ReorderBlocks(g2,JacobianIndices);
	  populateDerivativesContainer(g2,2,JacobianIndices);
	}

}
/* This version is not currently in use */
void KordpDynare::calcDerivatives(const Vector& yy, ogu::Jacobian& jacob)
{
//	ConstVector yym(yy, nstat(), nys());
//	ConstVector yyp(yy, nstat()+npred(), nyss());

//	Vector yyp(yy, nstat()+npred(), nyss());

	//double *g1, *g2;
	TwoDMatrix * jj= &jacob;
    Vector& out= *(new Vector(nY)); 
	out.zeros();
    Vector& xx= *(new Vector(nExog)); 
	xx.zeros();
	dynamicDLL.eval( yy,  xx, //int nb_row_x, 
				params, //int it_, 
				out, jj, NULL);
   //    model derivatives FSSparseTensor instance
        FSSparseTensor &mdTi=*(new FSSparseTensor (1, jj->ncols(),jj->nrows())); 
        for (int i = 0; i<jj->ncols(); i++){
                for (int j = 0; j<jj->nrows(); j++){
                    if (jj->get(i,j)!=0.0) // populate sparse if not zero
                        mdTi.insert(i, j, jj->get(i,j));
                }
        }
        // md container
        md.insert(&mdTi);
		delete &out;
		delete &xx;
}

void KordpDynare::calcDerivativesAtSteady()
{
	Vector xx(nexog());
	xx.zeros();
	calcDerivatives(*ySteady, xx);
}

/*******************************************************************************
* populateDerivatives to sparse Tensor and fit it in the Derivatives Container
*******************************************************************************/
void KordpDynare::populateDerivativesContainer(TwoDMatrix*g, int ord, const vector<int>*  vOrder)
{

#ifdef DEBUG
	mexPrintf("k_ord_dynare.cpp: populate FSSparseTensor in calcDerivatives: cols=%d , rows=%d\n"
        , g->ncols(),g->nrows());
#endif

	// model derivatives FSSparseTensor instance 
    FSSparseTensor *mdTi=(new FSSparseTensor (ord, nJcols,g->nrows()));

    IntSequence s(ord,0);
    s[0] = 0;
    s[1] = 0;
    int i = 0;
    while (i < g->ncols()){

      // insert new elements in each row
      if (ord == 1){
		for (int j = 0; j < g->nrows(); j++){
		  double x;
		  if (s[0] < nJcols-nExog)
			x = g->get(j,(*vOrder)[s[0]]);
		  else
			x = g->get(j,s[0]);
		  if (x != 0.0)
			mdTi->insert(s, j, x);
		}
		s[0]++;
      }
      else{	
		int s0, s1;
		if (s[0] < nJcols-nExog)
		  s0 = (*vOrder)[s[0]];
		else
		  s0 = s[0];
		if (s[1] < nJcols-nExog)
		  s1 = (*vOrder)[s[1]];
		else
		  s1 = s[1];
		if (s[1] >= s[0]){
		  s.print();
		  std::cout << s0 << " " << s1 << "\n";
		  int i1 = s0*nJcols+s1;
		  for (int j = 0; j < g->nrows(); j++){
			double x = g->get(j,i1);
			if (x != 0.0)
			  mdTi->insert(s, j, x);
		  }
		}
		s[1]++;
		// when one order index is finished
		// increase the previous one
		if (s[1] == nJcols){
		  s[0]++;
		  // update starting position of next indices
		  // in order to avoid dealing twice with the same
		  // symmetry. Increase matrix column counter
		  // accordingly
		  s[1] = 0;
		}
      }

      i++;

    }
    mdTi->print();
    // md container
    //md.clear();// this is to be used only for 1st order!!
    md.remove(Symmetry(ord));
    md.insert(mdTi);//(&mdTi);
}


void KordpDynare::writeModelInfo(Journal& jr) const
{
	// write info on variables
	{
		JournalRecordPair rp(journal);
		rp << "Information on variables" << endrec;
		JournalRecord rec1(journal);
		rec1 << "Number of endogenous:            " << ny() << endrec;
		JournalRecord rec2(journal);
		rec2 << "Number of exogenous:             " << nexog() << endrec;
		JournalRecord rec3(journal);
		rec3 << "Number of static:                " << nstat() << endrec;
		JournalRecord rec4(journal);
		rec4 << "Number of predetermined:         " << npred()+nboth() << endrec;
		JournalRecord rec5(journal);
		rec5 << "Number of forward looking:       " << nforw()+nboth() << endrec;
		JournalRecord rec6(journal);
		rec6 << "Number of both:                  " << nboth() << endrec;
	}

}
/*********************************************************
* LLxSteady()
* returns ySteady extended with leads and lags suitable for 
* passing to <model>_dynamic DLL
*************************************************************/
Vector * KordpDynare::LLxSteady( const Vector& yS){
	if ((nJcols-nExog) == yS.length()) {
		mexPrintf("k_ord_dynare.cpp: Warning in LLxSteady: ySteady already. right size");
		return NULL;
	}
	// create temporary square 2D matrix size nEndo x nEndo (sparse)  
	// for the lag, current and lead blocks of the jacobian
	Vector * llxSteady = new Vector(nJcols-nExog); 
	try{
		for (int ll_row=0; ll_row< ll_Incidence->nrows(); ll_row++)
		{
			// populate (non-sparse) vector with ysteady values 
			for (int i=0;i<nY;i++){
				if(ll_Incidence->get(ll_row,i))
					(*llxSteady)[((int)ll_Incidence->get(ll_row,i))-1] = yS[i];
			}
		}
	} catch (const TLException& e) {
		mexPrintf("Caugth TL exception in LLxSteady: ");
		e.print();
		return NULL;// 255;
	}catch (...){
		mexPrintf(" Error in LLxSteady - wrong index?");
	}

#ifdef DEBUG		
	for (int j=0;j<nJcols-nExog;j++)
			mexPrintf("LLxSteady: [%d] =%f .\n", j, (*llxSteady)[j]);
#endif	
	return llxSteady;
}


/************************************
* Reorder DynareJacobianIndices of variables in a vector according to
* given int * varOrder together with lead & lag incidence matrix and 
* any the extra columns for exogenous vars, and then, 
* reorders its blocks given by the varOrder and the Dynare++ expectations:

* extra	nboth+ npred (t-1) lags
* varOrder
		static:
		pred
		both    
		forward    
* extra both + nforw (t+1) leads, and 
* extra exogen

* so to match the jacobian organisation expected by the Appoximation class 
		both + nforw (t+1) leads
		static
		pred
		both
		forward
		nboth+ npred  (t-1) lags
		exogen
************************************/

vector<int> * KordpDynare::ReorderDynareJacobianIndices( const vector<int> * varOrder){
//	if ((nJcols != tdx->ncols()) && ( nY != tdx->nrows())) {
//		mexPrintf("k_ord_dynare.cpp: Error in ReorderBlocks: wrong size of jacobian");
//		return;
//	}
	// create temporary square 2D matrix size nEndo x nEndo (sparse)  
	// for the lag, current and lead blocks of the jacobian
//	int * JacobianIndices = (int*) calloc(nJcols+1, sizeof(int)); 
	vector<int> * JacobianIndices = new vector<int>(nJcols); 
	vector <int> tmp(nY); 
	int i,j, rjoff=nJcols-nExog-1; //, ll_off, j;
#ifdef DEBUG		
	mexPrintf("ReorderDynareJacobianIndice:ll_Incidence->nrows() =%d .\n", ll_Incidence->nrows());
#endif
	try{
		for (int ll_row=0; ll_row< ll_Incidence->nrows(); ll_row++)
		{
			// reorder in orde-var order & populate temporary nEndo (sparse) vector with 
			// the lag, current and lead blocks of the jacobian respectively
			for (i=0;i<nY;i++){
//				tmp[varOrder[j]]=(int)ll_Incidence->get(ll_row,j); 
				tmp[i]=((int)ll_Incidence->get(ll_row,(*varOrder)[i]-1)); 
#ifdef DEBUG		
	mexPrintf("get(ll_row,(*varOrder)[%d]-1)) = tmp[%d]=%d .\n", 
		i, i, (int)ll_Incidence->get(ll_row,(*varOrder)[i]-1));
#endif
			}
			// write the reordered blocks back to the jacobian
			// in reverse order
			for (j=nY-1;j>=0;j--){
				if (tmp[j]){
					(*JacobianIndices)[rjoff]=tmp[j] -1;
					rjoff--;
					if (rjoff<0){
	//					mexPrintf(" Error in ReorderIndices - negative rjoff index?");
						break;
	//					return NULL;
					}
				}
			}
		}
	} catch (const TLException& e) {
		mexPrintf("Caugth TL exception in ReorderIndices: ");
		e.print();
		return NULL;// 255;
	}catch (...){
		mexPrintf(" Error in ReorderIndices - wrong index?");
	}
	//add the indices for the nExog exogenous jacobians
	for (j=nJcols-nExog;j<nJcols;j++){
		(*JacobianIndices)[j]=j;
	}
#ifdef DEBUG		
	for (j=0;j<nJcols;j++)
			mexPrintf("ReorderDynareJacobianIndice: [%d] =%d .\n", j, (*JacobianIndices)[j]);
#endif	
	return JacobianIndices;
}


/************************************
* Reorder first set of columns of variables in a (jacobian) matrix 
* according to order given in  varsOrder together with the extras 
* assuming tdx ncols() - nExog is eaqual or less than length of varOrder and
* of any of its elements too.
************************************/

void KordpDynare::ReorderCols(TwoDMatrix * tdx, const vector<int> * vOrder){

	if (tdx->ncols() > vOrder->size()){
		mexPrintf(" Error in ReorderColumns - size of order var is too small");
		return;
	}
	TwoDMatrix tmp(*tdx); // temporary 2D matrix
	TwoDMatrix &tmpR=tmp;
	tdx->zeros();// empty original matrix
	// reorder the columns
	try{
		for (int i =0; i<tdx->ncols() ; i++)
			tdx->copyColumn(tmpR,(*vOrder)[i],i);
	} catch (const TLException& e) {
		printf("Caugth TL exception in ReorderColumns: ");
		e.print();
		return;// 255;
	}catch (...){
		mexPrintf(" Error in ReorderColumns - wrong index?");
	}
}
void KordpDynare::ReorderCols(TwoDMatrix * tdx, const int * vOrder){

	TwoDMatrix tmp(*tdx); // temporary 2D matrix
	TwoDMatrix &tmpR=tmp;
	tdx->zeros();// empty original matrix
	// reorder the columns
	try{
		for (int i =0; i<tdx->ncols() ; i++)
			tdx->copyColumn(tmpR,vOrder[i],i);
	} catch (const TLException& e) {
		printf("Caugth TL exception in ReorderColumns: ");
		e.print();
		return;// 255;
	}catch (...){
		mexPrintf(" Error in ReorderColumns - wrong index?");
	}
}

/***********************************************************************
* Recursive hierarchical block reordering of the higher order, input model 
*	derivatives inc. Hessian   
* This is now obsolete but kept in in case it is needed
***********************************************************************/

void KordpDynare::ReorderBlocks(TwoDMatrix * tdx, const vector<int> * vOrder){
	// determine order of the matrix

	double dbOrder = log(tdx->ncols())/log(nJcols);
	int ibOrder= (int) dbOrder;
	if ((double )ibOrder != dbOrder || ibOrder>nOrder) {
		mexPrintf(" Error in ReorderBlocks - wrong order %d", dbOrder);
		return;
	}
	TwoDMatrix tmp(*tdx); // temporary 2D matrix
	TwoDMatrix &tmpR=tmp;
	tdx->zeros();// empty original matrix

	if(ibOrder>1){
		int nBlocks=tmp.ncols()/ nJcols; //pow((float)nJcols,ibOrder-1);
		int bSize=tmp.ncols()/nBlocks;
		for (int j = 0; j<nBlocks;  ++j){
			TwoDMatrix subtdx(tmpR, bSize*((*vOrder)[j]), bSize);
			ReorderBlocks(&subtdx, vOrder);
			tdx->place(subtdx, 0, bSize*j);
		}
	} else{
		//ReorderColumns(TwoDMatrix * subtdx, const vector<int> * vOrder)
		if (tdx->ncols() > vOrder->size()){
			mexPrintf(" Error in ReorderColumns - size of order var is too small");
			return;
		}
		// reorder the columns
		try{
			for (int i =0; i<tdx->ncols() ; i++)
				tdx->copyColumn(tmpR,(*vOrder)[i],i);
		} catch (const TLException& e) {
			printf("Caugth TL exception in ReorderColumns: ");
			e.print();
			return;// 255;
		}catch (...){
			mexPrintf(" Error in ReorderColumns - wrong index?");
		}
	}
}


/****************************
*  K-Order Perturbation instance of Jacobian:
************************************/
KordpJacobian::KordpJacobian(KordpDynare& dyn)
	: Jacobian(dyn.ny()), dyn(dyn)
{
	zeros();
};

void KordpJacobian::eval(const Vector& yy)
{
		dyn.calcDerivatives( yy, *this);

};

void KordpVectorFunction::eval(const ConstVector& in, Vector& out)
{
	check_for_eval(in, out);
	Vector xx(d.nexog());
	xx.zeros();
	d.evaluateSystem(out, in, xx);
}

/**************************************************************************************/
/*       DynareNameList class                                                         */
/**************************************************************************************/
vector<int> DynareNameList::selectIndices(const vector<const char*>& ns) const
{
	vector<int> res;
	for (unsigned int i = 0; i < ns.size(); i++) {
		int j = 0;
		while (j < getNum() && strcmp(getName(j), ns[i]) != 0)
			j++;
		if (j == getNum())
			throw DynareException(__FILE__, __LINE__,
								  string("Couldn't find name for ") + ns[i] +
								  " in DynareNameList::selectIndices");
		res.push_back(j);
	}
	return res;
}

DynareNameList::DynareNameList(const  KordpDynare& dynare)
{
	for (int i = 0; i < dynare.ny(); i++) {
		names.push_back(dynare.dnl->getName(i));
	}
}
DynareNameList::DynareNameList(const KordpDynare& dynare, const char ** namesp)
{
#ifdef DEBUG		
    mexPrintf("k_ord_dynare DynareNameList.\n");
    mexPrintf("k_ord_dynare DynareNameList dynare.ny=%d .\n", dynare.ny());
#endif		
    
	for (int i = 0; i < dynare.ny(); i++) {
#ifdef DEBUG		
    mexPrintf("k_ord_dynare DynareNameList names[%d]= %s.\n", i, namesp[i] );
#endif		
		names.push_back(namesp[i]);
	}
}

DynareExogNameList::DynareExogNameList(const KordpDynare& dynare)
{
	for (int i = 0; i < dynare.nexog(); i++) {
		names.push_back(dynare.denl->getName(i));
	}
}

DynareExogNameList::DynareExogNameList(const KordpDynare& dynare, const char ** namesp)
{
#ifdef DEBUG		
    mexPrintf("k_ord_dynare DynareExogNameList dynare.nexog=%d .\n", dynare.nexog());
#endif		
	for (int i = 0; i < dynare.nexog(); i++) {
#ifdef DEBUG		
    mexPrintf("k_ord_dynare DynareExogNameList names[%d]= %s.\n", i, namesp[i] );
#endif		
		names.push_back(namesp[i]);
	}
}

DynareStateNameList::DynareStateNameList(const KordpDynare& dynare, const DynareNameList& dnl,
										 const DynareExogNameList& denl)
{
	for (int i = 0; i < dynare.nys(); i++){
#ifdef DEBUG		
    mexPrintf("k_ord_dynare DynareStateNameList dnl names[%d]= %s.\n", i, dnl.getName(i+dynare.nstat()) );
#endif		
		names.push_back(dnl.getName(i+dynare.nstat()));
    }
	for (int i = 0; i < dynare.nexog(); i++){
#ifdef DEBUG		
    mexPrintf("k_ord_dynare DynareStateNameList denl names[%d]= %s.\n", i, denl.getName(i));
#endif		
		names.push_back(denl.getName(i));
    }
}

