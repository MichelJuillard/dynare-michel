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
///#include "planner_builder.h"
///#include "forw_subst_builder.h"

#include "Dynare_pp/utils/cc/memory_file.h"
#include "Dynare_pp/utils/cc/exception.h"
///#include "parser/cc/parser_exception.h"
///#include "parser/cc/atom_substitutions.h"
#include "Dynare_pp/tl/cc/tl_exception.h"
#include "Dynare_pp/kord/kord_exception.h"

#ifndef DYNVERSION
#define DYNVERSION "unknown"
#endif


/**************************************************************************************/
/*       DynareNameList class                                                         */
/**************************************************************************************
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
*/////////////////////////
/**************************************************************************************/
/*       Dynare DynamicModel class                                                                 */
/**************************************************************************************/


KordpDynare::KordpDynare(const char** endo, int nstat,int npred, int nforw, int nboth
			   const char** exo, int nexog,
			   const char** par, int npar, Vector* ysteady,
			   const char* modName, int len, int order,
			   double sstol, Journal& jr)
	: journal(jr),  md(1)
{

	try{

	} else {
		throw DynareException(__FILE__, __LINE__, string("Could not open model file ")+modName);
	}

/// May need these later, GP, Oct. 08
///	dnl = new DynareNameList(*this);
///	denl = new DynareExogNameList(*this);
//	dsnl = new DynareStateNameList(*this, *dnl, *denl);
///	fe = new ogp::FormulaEvaluator(model->getParser());
///	fde = new ogp::FormulaDerEvaluator(model->getParser());
///	writeModelInfo(journal);
}

KordpDynare::KordpDynare(const KordpDynare& dynare)
	: ///journal(dynare.journal)///, model(NULL),
	  ysteady(NULL), md(dynare.md),
///	  dnl(NULL), denl(NULL), dsnl(NULL),
	ss_tol(dynare.ss_tol), nStat(dynare.nStat), nBoth(dynare.nBoth),
	nPred(dynare.nPred), nForw(dynare.nForw), nExo(dynare.nExo), 
	nYs(dynare.nYs), nYss(dynare.nYss),nY(dynare.nY), nOrder(dynare.nOrder)
{
///	model = dynare.model->clone();
	ysteady = new Vector(*(dynare.ysteady));
	params = new Vector(*(dynare.params));
	Vcov = new TwoDMatrix(*(dynare.Vcov));
	if (dynare.md)
		md = new TensorContainer<FSSparseTensor> (*(dynare.getModelDerivatives));
///	dnl = new DynareNameList(*this);
///	denl = new DynareExogNameList(*this);
///	dsnl = new DynareStateNameList(*this, *dnl, *denl);
}

KordpDynare::~KordpDynare()
{
	if (ysteady)
		delete ysteady;
//	if (params)
//		delete params;
//	if (Vcov)
//		delete Vcov;
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
	steady = (const Vector&) model->getInit();
	KordpDynareVectorFunction dvf(*this);
	KordpDynareJacobian dj(*this);
	ogu::NLSolver nls(dvf, dj, 500, ss_tol, journal);
	int iter;
	if (! nls.solve(steady, iter))
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

}

void KordpDynare::calcDerivatives(const Vector& yy, const Vector& xx)
{
	ConstVector yym(yy, nstat(), nys());
	ConstVector yyp(yy, nstat()+npred(), nyss());
	ogdyn::DynareAtomValues dav(model->getAtoms(), model->getParams(), yym, yy, yyp, xx);
	DynareDerEvalLoader ddel(model->getAtoms(), md, model->getOrder());
	for (int iord = 1; iord <= model->getOrder(); iord++)
		fde->eval(dav, ddel, iord);
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
}
/**********
DynareNameList::DynareNameList(const KordpDynare& dynare)
{
	for (int i = 0; i < dynare.ny(); i++) {
		int j = dynare.model->getAtoms().y2outer_endo()[i];
		const char* name = dynare.model->getAtoms().get_endovars()[j];
		names.push_back(name);
	}
}

DynareStateNameList::DynareStateNameList(const KordpDynare& dynare, const DynareNameList& dnl,
										 const DynareExogNameList& denl)
{
	for (int i = 0; i < dynare.nys(); i++)
		names.push_back(dnl.getName(i+dynare.nstat()));
	for (int i = 0; i < dynare.nexog(); i++)
		names.push_back(denl.getName(i));
}

DynareExogNameList::DynareExogNameList(const KordpDynare& dynare)
{
	for (int i = 0; i < dynare.nexog(); i++) {
		int j = dynare.model->getAtoms().y2outer_exo()[i];
		const char* name = dynare.model->getAtoms().get_exovars()[j];
		names.push_back(name);
	}
}
*****************///
/*****************
DynareEvalLoader::DynareEvalLoader(const ogp::FineAtoms& a, Vector& out)
	: Vector(out)
{
	if (a.ny() != out.length())
		throw DynareException(__FILE__, __LINE__, "Wrong length of out vector in DynareEvalLoader constructor");
}
**********/
/** This clears the container of model derivatives and initializes it
 * inserting empty sparse tensors up to the given order. */
 /*
DynareDerEvalLoader::DynareDerEvalLoader(const ogp::FineAtoms& a,
										 TensorContainer<FSSparseTensor>& mod_ders,
										 int order)
	: atoms(a), md(mod_ders)
{
	md.clear();
	for (int iord = 1; iord <= order; iord++) {
		FSSparseTensor* t = new FSSparseTensor(iord, atoms.ny()+atoms.nys()+atoms.nyss()+atoms.nexo(), atoms.ny());
		md.insert(t);
	}
}
void DynareDerEvalLoader::load(int i, int iord, const int* vars, double res)
{
	FSSparseTensor* t = md.get(Symmetry(iord));
	IntSequence s(iord, 0);
	for (int j = 0; j < iord; j++)
		s[j] = atoms.get_pos_of_all(vars[j]);
	t->insert(s, i, res);
}

DynareJacobian::DynareJacobian(Dynare& dyn)
	: Jacobian(dyn.ny()), d(dyn)
{
	zeros();
}

void DynareJacobian::eval(const Vector& yy)
{
	ogdyn::DynareSteadyAtomValues
		dav(d.getModel().getAtoms(), d.getModel().getParams(), yy);
	zeros();
	d.fde->eval(dav, *this, 1);
}

void DynareJacobian::load(int i, int iord, const int* vars, double res)
{
	if (iord != 1)
		throw DynareException(__FILE__, __LINE__,
							  "Derivative order different from order=1 in DynareJacobian::load");

	int t = vars[0];
	int j = d.getModel().getAtoms().get_pos_of_all(t);
	if (j < d.nyss())
		get(i, j+d.nstat()+d.npred()) += res;
	else if (j < d.nyss()+d.ny())
		get(i, j-d.nyss()) += res;
	else if (j < d.nyss()+d.ny()+d.nys())
		get(i, j-d.nyss()-d.ny()+d.nstat()) += res;
}
*/////////////////

/****************************
*  K-Order Perturbation instance of Jacobian:
************************************/
KordpJacobian::KordpJacobian(KordpDynare& dyn)
	: Jacobian(dyn.ny()), kdyn(dyn)
{
	zeros();
}

void KordpJacobian::eval(const Vector& yy)
{
///	ogdyn::DynareSteadyAtomValues
		dav(kdyn.getModel().getAtoms(), kdyn.getParams(), yy);
	zeros();
	d.fde->eval(dav, *this, 1);
}

void KordpJacobian::load(int i, int iord, const int* vars, double res)
{
	if (iord != 1)
		throw DynareException(__FILE__, __LINE__,
							  "Derivative order different from order=1 in KordpDynareJacobian::load");

	int t = vars[0];
	int j = d.getModel().getAtoms().get_pos_of_all(t);
	if (j < d.nyss())
		get(i, j+d.nstat()+d.npred()) += res;
	else if (j < d.nyss()+d.ny())
		get(i, j-d.nyss()) += res;
	else if (j < d.nyss()+d.ny()+d.nys())
		get(i, j-d.nyss()-d.ny()+d.nstat()) += res;
}
/////////////////


void KordpDynareVectorFunction::eval(const ConstVector& in, Vector& out)
{
	check_for_eval(in, out);
	Vector xx(d.nexog());
	xx.zeros();
	d.evaluateSystem(out, in, xx);
}

