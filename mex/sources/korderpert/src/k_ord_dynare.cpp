/* Dynare is free software: you can redistribute it and/or modify
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

#include "first_order.h"
#include "k_ord_dynare.h"

#include "mex.h" 

#include "memory_file.h"


//#include "k_order_perturbation.h"
#ifndef DYNVERSION
#define DYNVERSION "unknown"
#endif

    /******
void FistOrderApproximation::approxAtSteady()
{
	model.calcDerivativesAtSteady();
	FirstOrder fo(model.nstat(),model.npred(),model.nboth(),model.nforw(),
		model.nexog(),*(model.getModelDerivatives().get(Symmetry(1))),
		journal);
	KORD_RAISE_IF_X(!fo.isStable(),
		"The model is not Blanchard-Kahn stable",
		KORD_MD_NOT_STABLE);
	
	if(model.order()>=2){
		KOrder korder(model.nstat(),model.npred(),model.nboth(),model.nforw(),
			model.getModelDerivatives(),fo.getGy(),fo.getGu(),
			model.getVcov(),journal);
		korder.switchToFolded();
		for(int k= 2;k<=model.order();k++)
			korder.performStep<KOrder::fold> (k);
		
		saveRuleDerivs(korder.getFoldDers());
	}else{
		FirstOrderDerivs<KOrder::fold> fo_ders(fo);
		saveRuleDerivs(fo_ders);
        
	}
	check(0.0);
    Approximation::approxAtSteady();
    
	//saveRuleDerivs(fo);
}


void FistOrderApproximation::saveRuleDerivs(const FistOrder& fo)
{
	if(gy){
		delete gy;
		delete gu;
	}
	gy= new TwoDMatrix(fo.getGy);
	gu= new TwoDMatrix(fo.getGu);
}
     ****************/

/////////////////////////
/**************************************************************************************/
/*       Dynare DynamicModel class                                                                 */
/**************************************************************************************/
class KordpJacobian;

KordpDynare::KordpDynare(const char** endo,  int num_endo,
			   const char** exo, int nexog, int npar, //const char** par,
   			   Vector* ysteady, TwoDMatrix* vcov, Vector* inParams, int nstat,int npred, int nforw, int nboth,
			   const int nsteps, int norder, //const char* modName,
			   Journal& jr, DynamicModelDLL& dynamicDLL, double sstol)
	: nStat(nstat), nBoth(nboth), nPred(npred), nForw(nforw), nExog(nexog), nPar(npar),
	nYs(npred + nboth), nYss(nboth + nforw),nY(num_endo), nSteps(nsteps), nOrder(norder), journal(jr),  dynamicDLL(dynamicDLL),
	ySteady(ysteady), vCov(vcov), params (inParams),  md(1), dnl(NULL), denl(NULL), dsnl(NULL), ss_tol(sstol)
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

/****
		ySteady = new Vector(*(ySteady));
		params = new Vector(*(params));
		vCov = new TwoDMatrix(*(vCov));
******/
	//	throw DynareException(__FILE__, __LINE__, string("Could not open model file ")+modName);
	}
	catch (...){
        mexPrintf("k_ord_dynare: dynare constructor, error in StateNamelist construction.\n");
        throw DynareException(__FILE__, __LINE__, string("Could not construct Name Lists. \n"));
    }

/// May need these later, GP, Oct. 08
///	writeModelInfo(journal);
}

KordpDynare::KordpDynare(const KordpDynare& dynare)
	: nStat(dynare.nStat), nBoth(dynare.nBoth),	nPred(dynare.nPred), 
	nForw(dynare.nForw), nExog(dynare.nExog),  nPar(dynare.nPar),
	nYs(dynare.nYs), nYss(dynare.nYss),nY(dynare.nY), 
	nSteps(dynare.nSteps), nOrder(dynare.nOrder), journal(dynare.journal),
	dynamicDLL(dynare.dynamicDLL), //modName(dynare.modName),
	ySteady(NULL), params(NULL), vCov(NULL), md(dynare.md), 
	dnl(NULL), denl(NULL), dsnl(NULL), ss_tol(dynare.ss_tol)
{
///	model = dynare.model->clone();
	ySteady = new Vector(*(dynare.ySteady));
	params = new Vector(*(dynare.params));
	vCov = new TwoDMatrix(*(dynare.vCov));
//	if (dynare.md)// !=NULL)
//	Inititalise ModelDerivatives md
//	md= *(new TensorContainer<FSSparseTensor>(dynare.md));
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
/***
	ConstVector yym(yy, nstat(), nys());
	ConstVector yyp(yy, nstat()+npred(), nyss());
	evaluateSystem(out, yym, yy, yyp, xx);
***/
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
/*//////////////////////////
	ogdyn::DynareAtomValues dav(model->getAtoms(), model->getParams(), yym, yy, yyp, xx);
	DynareEvalLoader del(model->getAtoms(), out);
	fe->eval(dav, del);
///////////////////////*/
#ifdef DEBUG
	mexPrintf("k_order_dynaare.cpp: Call eval in EvaluateSystem\n");
#endif
	dynamicDLL.eval( yy,  xx, //int nb_row_x, 
				 params, //int it_, 
                 out, NULL, NULL);
}
//void KordpDynare::calcDerivatives(const Vector& yy, const TwoDMatrix& xx)
void KordpDynare::calcDerivatives(const Vector& yy, const Vector& xx)
{
//	ConstVector yym(yy, nstat(), nys());
//	ConstVector yyp(yy, nstat()+npred(), nyss());

//	Vector yyp(yy, nstat()+npred(), nyss());

	//double *g1, *g2;
    TwoDMatrix *g1;//, *g2;
	g1=new TwoDMatrix(0,0); // just a signal: generate something so g1 (and out) are not null
    Vector& out= *(new Vector(nY));
#ifdef DEBUG
	mexPrintf("k_order_dynaare.cpp: Call eval in calcDerivatives\n");
#endif
	dynamicDLL.eval( yy,  xx, //int nb_row_x, 
				params, //int it_, 
				out, g1, NULL);
#ifdef DEBUG
	mexPrintf("k_order_dynaare.cpp: populate FSSparseTensor in calcDerivatives: cols=%d , rows=%d\n"
        , g1->ncols(),g1->nrows());
#endif

   //    model derivatives FSSparseTensor instance for single order only 
    //(higher orders requires Symetry to insert in particular position.)
        FSSparseTensor mdTi=*(new FSSparseTensor (1, g1->ncols(),g1->nrows())); 
        for (int i = 0; i<g1->ncols(); i++){
                for (int j = 0; j<g1->nrows(); j++){
                    if (g1->get(i,j)!=0.0) // populate sparse if not zero
                        mdTi.insert(i, j,g1->get(i,j));
                }
        }
        // md container
//        md=*(new TensorContainer<FSSparseTensor>(1)); 
        md.clear();
        md.insert(&mdTi);
}


void KordpDynare::calcDerivatives(const Vector& yy, ogu::Jacobian& jacob)
{
//	ConstVector yym(yy, nstat(), nys());
//	ConstVector yyp(yy, nstat()+npred(), nyss());

//	Vector yyp(yy, nstat()+npred(), nyss());

	//double *g1, *g2;
	TwoDMatrix * jj= &jacob;
    Vector& out= *(new Vector(nY));
    Vector& xx= *(new Vector(nExog));
	dynamicDLL.eval( yy,  xx, //int nb_row_x, 
				params, //int it_, 
				out, jj, NULL);
   //    model derivatives FSSparseTensor instance
        FSSparseTensor mdTi=*(new FSSparseTensor (1, jj->ncols(),jj->nrows())); 
        for (int i = 0; i<jj->ncols(); i++){
                for (int j = 0; j<jj->nrows(); j++){
                    if (jj->get(i,j)!=0.0) // populate sparse if not zero
                        mdTi.insert(i, j, jj->get(i,j));
                }
        }
        // md container
//        md=*(new TensorContainer<FSSparseTensor>(1)); 
        md.clear();
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

/************** May be needed
	// write info on forward substitutions
	const ogdyn::ForwSubstInfo* finfo = model->get_forw_subst_info();
	if (finfo) {
		JournalRecordPair rp(journal);
		rp << "Information on forward substitutions" << endrec;
		JournalRecord rec1(journal);
		rec1 << "Number of affected equations:    " << finfo->num_affected_equations << endrec;
		JournalRecord rec2(journal);
		rec2 << "Number of substituted terms:     " << finfo->num_subst_terms << endrec;
		JournalRecord rec3(journal);
		rec3 << "Number of auxiliary variables:   " << finfo->num_aux_variables << endrec;
		JournalRecord rec4(journal);
		rec4 << "Number of new terms in the tree: " << finfo->num_new_terms << endrec;
	}

	// write info on substitutions
	const ogp::SubstInfo* sinfo = model->get_subst_info();
	if (sinfo) {
		JournalRecordPair rp(journal);
		rp << "Information on substitutions" << endrec;
		JournalRecord rec1(journal);
		rec1 << "Number of substitutions:         " << sinfo->num_substs << endrec;
	}
	**************/
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

