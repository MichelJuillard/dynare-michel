
#include "dynare3.h"
#include "dynare_exception.h"
#include "planner_builder.h"
#include "forw_subst_builder.h"

#include "utils/cc/memory_file.h"
#include "utils/cc/exception.h"
#include "parser/cc/parser_exception.h"
#include "parser/cc/atom_substitutions.h"
#include "../tl/cc/tl_exception.h"
#include "../kord/kord_exception.h"

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

/**************************************************************************************/
/*       Dynare class                                                                 */
/**************************************************************************************/

Dynare::Dynare(const char* modname, int ord, double sstol, Journal& jr)
	: journal(jr), model(NULL), ysteady(NULL), md(1), dnl(NULL), denl(NULL), dsnl(NULL),
	  fe(NULL), fde(NULL), ss_tol(sstol)
{
	// make memory file
	ogu::MemoryFile mf(modname);
	if (mf.exists()) {
		try {
			model = new ogdyn::DynareParser(mf.base(), mf.length(), ord);
		} catch (const ogp::ParserException& pe) {
			int line;
			int col;
			mf.line_and_col(pe.offset(), line, col);
			throw DynareException(pe.message(), modname, line, col);
		}
		ysteady = new Vector(model->getAtoms().ny());
		dnl = new DynareNameList(*this);
		denl = new DynareExogNameList(*this);
		dsnl = new DynareStateNameList(*this, *dnl, *denl);
		fe = new ogp::FormulaEvaluator(model->getParser());
		fde = new ogp::FormulaDerEvaluator(model->getParser());
		writeModelInfo(journal);
	} else {
		throw DynareException(__FILE__, __LINE__, string("Could not open model file ")+modname);
	}
}

Dynare::Dynare(const char** endo, int num_endo,
			   const char** exo, int num_exo,
			   const char** par, int num_par,
			   const char* equations, int len, int ord,
			   double sstol, Journal& jr)
	: journal(jr), model(NULL), ysteady(NULL), md(1), dnl(NULL), denl(NULL), dsnl(NULL),
	  fe(NULL), fde(NULL), ss_tol(sstol)
{
	try {
		model = new ogdyn::DynareSPModel(endo, num_endo, exo, num_exo, par, num_par,
										 equations, len, ord);
	} catch (const ogp::ParserException& pe) {
		throw DynareException(pe.message(), pe.offset());
	}
	ysteady = new Vector(model->getAtoms().ny());
	dnl = new DynareNameList(*this);
	denl = new DynareExogNameList(*this);
	dsnl = new DynareStateNameList(*this, *dnl, *denl);
	fe = new ogp::FormulaEvaluator(model->getParser());
	fde = new ogp::FormulaDerEvaluator(model->getParser());
	writeModelInfo(journal);
}

Dynare::Dynare(const Dynare& dynare)
	: journal(dynare.journal), model(NULL),
	  ysteady(NULL), md(dynare.md),
	  dnl(NULL), denl(NULL), dsnl(NULL), fe(NULL), fde(NULL),
	  ss_tol(dynare.ss_tol)
{
	model = dynare.model->clone();
	ysteady = new Vector(*(dynare.ysteady));
	dnl = new DynareNameList(*this);
	denl = new DynareExogNameList(*this);
	dsnl = new DynareStateNameList(*this, *dnl, *denl);
	fe = new ogp::FormulaEvaluator(model->getParser());
	fde = new ogp::FormulaDerEvaluator(model->getParser());
}

Dynare::~Dynare()
{
	if (model)
		delete model;
	if (ysteady)
		delete ysteady;
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
}

void Dynare::writeMat(mat_t* fd, const char* prefix) const
{
	char tmp[100];
	sprintf(tmp, "%s_vars", prefix);
	getAllEndoNames().writeMat(fd, tmp);
	getAllEndoNames().writeMatIndices(fd, prefix);
	sprintf(tmp, "%s_state_vars", prefix);
	getStateNames().writeMat(fd, tmp);
	sprintf(tmp, "%s_shocks", prefix);
	getExogNames().writeMat(fd, tmp);
	getExogNames().writeMatIndices(fd, prefix);
	sprintf(tmp, "%s_vcov_exo", prefix);
	model->getVcov().writeMat(fd, tmp);
	TwoDMatrix aux(1,1);
	sprintf(tmp, "%s_nstat", prefix);
	aux.get(0,0) = nstat();
	aux.writeMat(fd, tmp);
	sprintf(tmp, "%s_npred", prefix);
	aux.get(0,0) = npred();
	aux.writeMat(fd, tmp);
	sprintf(tmp, "%s_nboth", prefix);
	aux.get(0,0) = nboth();
	aux.writeMat(fd, tmp);
	sprintf(tmp, "%s_nforw", prefix);
	aux.get(0,0) = nforw();
	aux.writeMat(fd, tmp);
}

void Dynare::writeDump(const std::string&  basename) const
{
	std::string fname(basename);
	fname += ".dump";
	std::ofstream out(fname.c_str());
	model->dump_model(out);
	out.close();
}

void Dynare::solveDeterministicSteady(Vector& steady)
{
	JournalRecordPair pa(journal);
	pa << "Non-linear solver for deterministic steady state" << endrec;
	steady = (const Vector&) model->getInit();
	DynareVectorFunction dvf(*this);
	DynareJacobian dj(*this);
	ogu::NLSolver nls(dvf, dj, 500, ss_tol, journal);
	int iter;
	if (! nls.solve(steady, iter))
		throw DynareException(__FILE__, __LINE__,
							  "Could not obtain convergence in non-linear solver");
}

// evaluate system at given y_t=y_{t+1}=y_{t-1}, and given shocks x_t
void Dynare::evaluateSystem(Vector& out, const Vector& yy, const Vector& xx)
{
	ConstVector yym(yy, nstat(), nys());
	ConstVector yyp(yy, nstat()+npred(), nyss());
	evaluateSystem(out, yym, yy, yyp, xx);
}

// evaluate system at given y^*_{t-1}, y_t, y^{**}_{t+1} and at
// exogenous x_t, all three vectors yym, yy, and yyp have the
// respective lengths of y^*_{t-1}, y_t, y^{**}_{t+1}
void Dynare::evaluateSystem(Vector& out, const Vector& yym, const Vector& yy,
							const Vector& yyp, const Vector& xx)
{
	ogdyn::DynareAtomValues dav(model->getAtoms(), model->getParams(), yym, yy, yyp, xx);
	DynareEvalLoader del(model->getAtoms(), out);
	fe->eval(dav, del);
}

void Dynare::calcDerivatives(const Vector& yy, const Vector& xx)
{
	ConstVector yym(yy, nstat(), nys());
	ConstVector yyp(yy, nstat()+npred(), nyss());
	ogdyn::DynareAtomValues dav(model->getAtoms(), model->getParams(), yym, yy, yyp, xx);
	DynareDerEvalLoader ddel(model->getAtoms(), md, model->getOrder());
	for (int iord = 1; iord <= model->getOrder(); iord++)
		fde->eval(dav, ddel, iord);
}

void Dynare::calcDerivativesAtSteady()
{
	Vector xx(nexog());
	xx.zeros();
	calcDerivatives(*ysteady, xx);
}

void Dynare::writeModelInfo(Journal& jr) const
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

	// write info on planner variables
	const ogdyn::PlannerInfo* pinfo = model->get_planner_info();
	if (pinfo) {
		JournalRecordPair rp(journal);
		rp << "Information on planner variables" << endrec;
		JournalRecord rec1(journal);
		rec1 << "Number of Lagrange multipliers:  " << pinfo->num_lagrange_mults << endrec;
		JournalRecord rec2(journal);
		rec2 << "Number of auxiliary variables:   " << pinfo->num_aux_variables << endrec;
		JournalRecord rec3(journal);
		rec3 << "Number of new terms in the tree: " << pinfo->num_new_terms << endrec;
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

DynareNameList::DynareNameList(const Dynare& dynare)
{
	for (int i = 0; i < dynare.ny(); i++) {
		int j = dynare.model->getAtoms().y2outer_endo()[i];
		const char* name = dynare.model->getAtoms().get_endovars()[j];
		names.push_back(name);
	}
}

DynareStateNameList::DynareStateNameList(const Dynare& dynare, const DynareNameList& dnl,
										 const DynareExogNameList& denl)
{
	for (int i = 0; i < dynare.nys(); i++)
		names.push_back(dnl.getName(i+dynare.nstat()));
	for (int i = 0; i < dynare.nexog(); i++)
		names.push_back(denl.getName(i));
}

DynareExogNameList::DynareExogNameList(const Dynare& dynare)
{
	for (int i = 0; i < dynare.nexog(); i++) {
		int j = dynare.model->getAtoms().y2outer_exo()[i];
		const char* name = dynare.model->getAtoms().get_exovars()[j];
		names.push_back(name);
	}
}

DynareEvalLoader::DynareEvalLoader(const ogp::FineAtoms& a, Vector& out)
	: Vector(out)
{
	if (a.ny() != out.length())
		throw DynareException(__FILE__, __LINE__, "Wrong length of out vector in DynareEvalLoader constructor");
}

/** This clears the container of model derivatives and initializes it
 * inserting empty sparse tensors up to the given order. */
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

void DynareVectorFunction::eval(const ConstVector& in, Vector& out)
{
	check_for_eval(in, out);
	Vector xx(d.nexog());
	xx.zeros();
	d.evaluateSystem(out, in, xx);
}

