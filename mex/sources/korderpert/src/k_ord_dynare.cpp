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
// based on: work by O.Kamenik


#include "k_ord_dynare.h"
#include "dynare_exception.h"
#include "decision_rule.h"
#include "fs_tensor.h"
#include "SylvException.h"

#include "mex.h" 

///#include "forw_subst_builder.h"

#include "Dynare_pp/utils/cc/memory_file.h"
#include "Dynare_pp/utils/cc/exception.h"
#include "Dynare_pp/tl/cc/tl_exception.h"
#include "Dynare_pp/kord/kord_exception.h"


//#include "k_order_perturbation.h"
#ifndef DYNVERSION
#define DYNVERSION "unknown"
#endif


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
/////////////////////////
/**************************************************************************************/
/*       Dynare DynamicModel class                                                                 */
/**************************************************************************************/


KordpDynare::KordpDynare(const char** endo,  int num_endo,
			   const char** exo, int nexog,
			   const char** par, int nPar, 
   			   Vector* ysteady, int nstat,int nPred, int nForw, int nboth,
			   const char* modName, int len, int nOrder,
			   double sstol, Journal& jr, DynamicModelDLL& dynamicDLL)
	: journal(jr),  md(1), 	  ysteady(NULL),
	  dnl(NULL), denl(NULL), dsnl(NULL), nPar(nPar), nOrder(nOrder),
	ss_tol(sstol), nStat(nstat), nBoth(nboth),
	 nForw(nForw), nExog(nexog), nPred(nPred),
	nYs(nYs), nYss(nYss),nY(nY), dynamicDLL(dynamicDLL)

{

	try{

		dnl = new DynareNameList(*this);
		denl = new DynareExogNameList(*this);
		dsnl = new DynareStateNameList(*this, *dnl, *denl);

		throw DynareException(__FILE__, __LINE__, string("Could not open model file ")+modName);
	}
	catch (...)
	{}

/// May need these later, GP, Oct. 08
///	writeModelInfo(journal);
}

class KordpJacobian;

KordpDynare::KordpDynare(const KordpDynare& dynare)
	: journal(dynare.journal),///, model(NULL),
	  ysteady(NULL), md(dynare.md),
	  dnl(NULL), denl(NULL), dsnl(NULL),
	ss_tol(dynare.ss_tol), nStat(dynare.nStat), nBoth(dynare.nBoth),
	nPred(dynare.nPred), nForw(dynare.nForw), nExog(dynare.nExog), 
	nYs(dynare.nYs), nYss(dynare.nYss),nY(dynare.nY),
	dynamicDLL(dynare.dynamicDLL),nOrder(dynare.nOrder), nPar(dynare.nPar)
{
///	model = dynare.model->clone();
	ysteady = new Vector(*(dynare.ysteady));
	params = new Vector(*(dynare.params));
	Vcov = new TwoDMatrix(*(dynare.Vcov));
//	if (dynare.md)// !=NULL)
//	Inititalise ModelDerivatives md
	md= *(new TensorContainer<FSSparseTensor>::TensorContainer(dynare.md));
	dnl = new DynareNameList(*this);
	denl = new DynareExogNameList(*this);
	dsnl = new DynareStateNameList(*this, *dnl, *denl);
}

KordpDynare::~KordpDynare()
{
	if (ysteady)
		delete ysteady;
	if (params)
		delete params;
	if (Vcov)
		delete Vcov;
/*************** May be needed
	if (dnl)
		delete dnl;
	if (dsnl)
		delete dsnl;
	if (denl)
		delete denl;
	if (fe)
		delete fe;
	if (fde)
		delete fde;
***********///

}

void KordpDynare::solveDeterministicSteady(Vector& steady)
{
	JournalRecordPair pa(journal);
	pa << "Non-linear solver for deterministic steady state" << endrec;
	//steady = (const Vector&) model->getInit();
	KordpVectorFunction dvf(*this);
	KordpJacobian dj(*this);
	ogu::NLSolver nls(dvf, dj, 500, ss_tol, journal);
	int iter;
	if (! nls.solve(*ysteady, iter))
		throw DynareException(__FILE__, __LINE__,
							  "Could not obtain convergence in non-linear solver");
}

// evaluate system at given y_t=y_{t+1}=y_{t-1}, and given shocks x_t
void KordpDynare::evaluateSystem(Vector& out, const Vector& yy, const Vector& xx)
{
	ConstVector yym(yy, nstat(), nys());
	ConstVector yyp(yy, nstat()+npred(), nyss());
	evaluateSystem(out, yym, yy, yyp, xx);
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
	mexPrintf("Call in EvaluateSystem\n");
#endif
//	DynamicDLL->eval(double *y, double *x, int nb_row_x, double *params, int it_, 
//				double *residual, double *g1, double *g2);
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
    TwoDMatrix *g1, *g2;
    Vector& out= *(new Vector(ny()));
	dynamicDLL.eval( yy,  xx, //int nb_row_x, 
				params, //int it_, 
				out, g1, NULL);

//    if (!md){
//        md=new TensorContainer<FSSparseTensor> (1); 
    //    model derivatives FSSparseTensor instance
        FSSparseTensor mdTi=*(new FSSparseTensor (1, g1->ncols(),g1->nrows())); 
        for (int i = 0; i<g1->ncols(); i++){
                for (int j = 0; j<g1->nrows(); j++){
                    if (g1->get(i,j)!=0.0) // populate sparse if not zero
                        mdTi.insert(i, j,g1->get(i,j));
                }
        }
        // md container
        md=*(new TensorContainer<FSSparseTensor>::TensorContainer(1)); 

        md.insert(&mdTi);
//    }else{}
		
}
void KordpDynare::calcDerivativesAtSteady()
{
	Vector xx(nexog());
	xx.zeros();
	calcDerivatives(*ysteady, xx);
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

