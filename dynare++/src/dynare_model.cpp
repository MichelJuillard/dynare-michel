// Copyright (C) 2006, Ondra Kamenik

// $Id: dynare_model.cpp 2269 2008-11-23 14:33:22Z michel $

#include "parser/cc/parser_exception.h"
#include "parser/cc/location.h"
#include "utils/cc/exception.h"
#include "dynare_model.h"
#include "dynare_exception.h"
#include "planner_builder.h"
#include "forw_subst_builder.h"

#include <cstdlib>

#include <string>
#include <cmath>
#include <climits>

using namespace ogdyn;

ParsedMatrix::ParsedMatrix(const ogp::MatrixParser& mp)
	: TwoDMatrix(mp.nrows(), mp.ncols())
{
	zeros();
	for (ogp::MPIterator it = mp.begin(); it != mp.end(); ++it)
		get(it.row(), it.col()) = *it;
}

DynareModel::DynareModel()
	: atoms(), eqs(atoms), order(-1),
	  param_vals(0), init_vals(0), vcov_mat(0),
	  t_plobjective(-1), t_pldiscount(-1),
	  pbuilder(NULL), fbuilder(NULL),
	  atom_substs(NULL), old_atoms(NULL)
{}

DynareModel::DynareModel(const DynareModel& dm)
	: atoms(dm.atoms), eqs(dm.eqs, atoms), order(dm.order),
	  param_vals(0), init_vals(0), vcov_mat(0),
	  t_plobjective(dm.t_plobjective),
	  t_pldiscount(dm.t_pldiscount),
	  pbuilder(NULL), fbuilder(NULL),
	  atom_substs(NULL), old_atoms(NULL)
{
	if (dm.param_vals)
		param_vals = new Vector((const Vector&)*(dm.param_vals));
	if (dm.init_vals)
		init_vals = new Vector((const Vector&)*(dm.init_vals));
	if (dm.vcov_mat)
		vcov_mat = new TwoDMatrix((const TwoDMatrix&)*(dm.vcov_mat));
	if (dm.old_atoms)
		old_atoms = new DynareDynamicAtoms((const DynareDynamicAtoms&)*(dm.old_atoms));
	if (dm.atom_substs)
		atom_substs = new ogp::AtomSubstitutions((*dm.atom_substs), *old_atoms, atoms);
	if (dm.pbuilder)
		pbuilder = new PlannerBuilder(*(dm.pbuilder), *this);
	if (dm.fbuilder)
		fbuilder = new ForwSubstBuilder(*(dm.fbuilder), *this);
}

DynareModel::~DynareModel()
{
	if (param_vals)
		delete param_vals;
	if (init_vals)
		delete init_vals;
	if (vcov_mat)
		delete vcov_mat;
	if (old_atoms)
		delete old_atoms;
	if (atom_substs)
		delete atom_substs;
	if (pbuilder)
		delete pbuilder;
	if (fbuilder)
		delete fbuilder;
}

const PlannerInfo* DynareModel::get_planner_info() const
{
	if (pbuilder)
		return &(pbuilder->get_info());
	return NULL;
}

const ForwSubstInfo* DynareModel::get_forw_subst_info() const
{
	if (fbuilder)
		return &(fbuilder->get_info());
	return NULL;
}

const ogp::SubstInfo* DynareModel::get_subst_info() const
{
	if (atom_substs)
		return &(atom_substs->get_info());
	return NULL;
}

void DynareModel::setInitOuter(const Vector& x)
{
	if (x.length() != atoms.ny())
		throw DynareException(__FILE__, __LINE__,
							  "Wrong length of vector in DynareModel::setInitOuter");
	for (int i = 0; i < atoms.ny(); i++)
		(*init_vals)[i] = x[atoms.y2outer_endo()[i]];
}

void DynareModel::print() const
{
	printf("all atoms:\n");
	atoms.print();
	printf("formulas:\n");
	DebugOperationFormatter dof(*this);
	for (int i = 0; i < eqs.nformulas(); i++) {
		int tf = eqs.formula(i);
		printf("formula %d:\n", tf);
		eqs.getTree().print_operation_tree(tf, stdout, dof);
	}
}

void DynareModel::dump_model(std::ostream& os) const
{
	// endogenous variable declaration
	os << "var";
	for (int i = 0; i < (int)atoms.get_endovars().size(); i++)
		os << " " << atoms.get_endovars()[i];
	os << ";\n\n";

	// exogenous variables
	os << "varexo";
	for (int i = 0; i < (int)atoms.get_exovars().size(); i++)
		os << " " << atoms.get_exovars()[i];
	os << ";\n\n";

	// parameters
	os << "parameters";
	for (int i = 0; i < (int)atoms.get_params().size(); i++)
		os << " " << atoms.get_params()[i];
	os << ";\n\n";

	// parameter values
	os.precision(16);
	for (int i = 0; i < (int)atoms.get_params().size(); i++)
		os << atoms.get_params()[i] << "=" << getParams()[i] << ";\n";
	os << "\n\n";

	// model section
	ogp::OperationStringConvertor osc(atoms, getParser().getTree());
	os << "model;\n";
	for (int i = 0; i < getParser().nformulas(); i++) {
		os << "// Equation " << i << "\n0 = ";
		int t = getParser().formula(i);
		os << osc.convert(getParser().getTree().operation(t), t);
		os << ";\n";
	}
	os << "end;\n";

	// initval as steady state
	os << "initval;\n";
	for (int i = 0; i < (int)atoms.get_endovars().size(); i++)
		os << atoms.get_endovars()[atoms.y2outer_endo()[i]] << "=" << getInit()[i] << ";\n";
	os << "end;\n";
}

void DynareModel::add_name(const char* name, int flag)
{
	if (flag == 1) {
		// endogenous
		atoms.register_uniq_endo(name);
	} else if (flag == 2) {
		// exogenous
		atoms.register_uniq_exo(name);
	} else if (flag == 3) {
		// parameter
		atoms.register_uniq_param(name);
	} else {
		throw DynareException(__FILE__, __LINE__,
							  "Unrecognized flag value.");
	}
}

void DynareModel::check_model() const
{
	if (order == -1)
		throw DynareException(__FILE__,__LINE__,
							  "Order of approximation not set in DynareModel::check_model");

	if (atoms.ny() != eqs.nformulas()) {
		char mes[1000];
		sprintf(mes, "Model has %d equations for %d endogenous variables", eqs.nformulas(), atoms.ny());
		throw DynareException(__FILE__, __LINE__, mes);
	}
	
	// check whether all nulary terms of all formulas in eqs are
	// either constant or assigned to a name
	for (int i = 0; i < eqs.nformulas(); i++) {
		int ft = eqs.formula(i);
		const hash_set<int>& nuls = eqs.nulary_of_term(ft);
		for (hash_set<int>::const_iterator it = nuls.begin();
			 it != nuls.end(); ++it)
			if (! atoms.is_constant(*it) && ! atoms.is_named_atom(*it))
				throw DynareException(__FILE__,__LINE__,
									  "Dangling nulary term found, internal error.");
	}

	int mlag, mlead;
	atoms.exovarspan(mlead, mlag);
	if (atoms.nexo() > 0 && (mlead != 0 || mlag != 0))
		throw DynareException(__FILE__,__LINE__,
							  "The model contains occurrences of lagged/leaded exogenous variables");

	atoms.endovarspan(mlead, mlag);
	if (mlead > 1 || mlag < -1)
		throw DynareException(__FILE__,__LINE__,
							  "The model contains occurrences of too lagged/leaded endogenous variables");

	// check the dimension of vcov matrix
	if (getAtoms().nexo() != getVcov().nrows())
		throw DynareException(__FILE__,__LINE__,
							  "Dimension of VCOV matrix does not correspond to the shocks");
}

int DynareModel::variable_shift(int t, int tshift)
{
	const char* name = atoms.name(t);
	if (atoms.is_type(name, DynareDynamicAtoms::param) ||
		atoms.is_constant(t))
		throw DynareException(__FILE__, __LINE__,
							  "The tree index is not a variable in DynareModel::variable_shift");
	int ll = atoms.lead(t) + tshift;
	int res = atoms.index(name, ll);
	if (res == -1) {
		std::string str(name);
		str += '(';
		char tmp[50];
		sprintf(tmp,"%d",ll);
		str += tmp;
		str += ')';
		res = eqs.add_nulary(str.c_str());
	}
	return res;
}

void DynareModel::variable_shift_map(const hash_set<int>& a_set, int tshift,
									 map<int,int>& s_map)
{
	s_map.clear();
	for (hash_set<int>::const_iterator it = a_set.begin();
		 it != a_set.end(); ++it) {
		int t = *it;
		// make shift map only for non-constants and non-parameters
		if (! atoms.is_constant(t)) {
			const char* name = atoms.name(t);
			if (atoms.is_type(name, DynareDynamicAtoms::endovar) ||
				atoms.is_type(name, DynareDynamicAtoms::exovar)) {
				int tt = variable_shift(t, tshift);
				s_map.insert(map<int,int>::value_type(t,tt));
			}
		}
	}
}

void DynareModel::termspan(int t, int& mlead, int& mlag) const
{
	mlead = INT_MIN;
	mlag = INT_MAX;
	const hash_set<int>& nul_terms = eqs.nulary_of_term(t);
	for (hash_set<int>::const_iterator ni = nul_terms.begin();
		 ni != nul_terms.end(); ++ni) {
		if (!atoms.is_constant(*ni) &&
			(atoms.is_type(atoms.name(*ni), DynareDynamicAtoms::endovar) ||
			 atoms.is_type(atoms.name(*ni), DynareDynamicAtoms::exovar))) {
			int ll = atoms.lead(*ni);
			if (ll < mlag)
				mlag = ll;
			if (ll > mlead)
				mlead = ll;
		}
	}
}

bool DynareModel::is_constant_term(int t) const
{
	const hash_set<int>& nul_terms = eqs.nulary_of_term(t);
	for (hash_set<int>::const_iterator ni = nul_terms.begin();
		 ni != nul_terms.end(); ++ni)
		if (! atoms.is_constant(*ni) &&
			! atoms.is_type(atoms.name(*ni), DynareDynamicAtoms::param))
			return false;
	return true;
}

hash_set<int> DynareModel::get_nonlinear_subterms(int t) const
{
	NLSelector nls(*this);
	return eqs.getTree().select_terms(t, nls);
}

void DynareModel::substitute_atom_for_term(const char* name, int ll, int t)
{
	// if the term t is itself a named atom (parameter, exo, endo),
	// then we have to unassign it first
	if (atoms.is_named_atom(t))
		atoms.unassign_variable(atoms.name(t), atoms.lead(t), t);
	// assign allocated tree index
	// for the term now to name(ll)
	atoms.assign_variable(name, ll, t);
	// make operation t nulary in operation tree
	eqs.nularify(t);
}

void DynareModel::final_job()
{
	if (t_plobjective != -1 && t_pldiscount != -1) {
		// at this moment include all equations and all variables; in
		// future we will exclude purely exogenous processes; todo:
		PlannerBuilder::Tvarset vset;
		for (int i = 0; i < atoms.ny(); i++)
			vset.insert(atoms.get_endovars()[i]);
		PlannerBuilder::Teqset eset;
		for (int i = 0; i < eqs.nformulas(); i++)
			eset.push_back(i);

		// construct the planner builder, this adds a lot of stuff to
		// the model
		if (pbuilder)
			delete pbuilder;
		pbuilder = new PlannerBuilder(*this, vset, eset);
	}

	// construct ForwSubstBuilder
	if (fbuilder)
		delete fbuilder;
	fbuilder = new ForwSubstBuilder(*this);

	// call parsing_finished (this will define an outer ordering of all variables)
	atoms.parsing_finished(ogp::VarOrdering::bfspbfpb);
    // make a copy of atoms and name it old_atoms
	if (old_atoms)
		delete old_atoms;
	old_atoms = new DynareDynamicAtoms(atoms);
	// construct empty substitutions from old_atoms to atoms
	if (atom_substs)
		delete atom_substs;
	atom_substs = new ogp::AtomSubstitutions(*old_atoms, atoms);
	// do the actual substitution, it will also call
	// parsing_finished for atoms which creates internal orderings
	atoms.substituteAllLagsAndExo1Leads(eqs, *atom_substs);
}

extern ogp::location_type dynglob_lloc;

DynareParser::DynareParser(const char* stream, int len, int ord)
	: DynareModel(),
	  pa_atoms(), paramset(pa_atoms),
      ia_atoms(), initval(ia_atoms), vcov(),
	  model_beg(0), model_end(-1),
	  paramset_beg(0), paramset_end(-1),
	  initval_beg(0), initval_end(-1),
	  vcov_beg(0), vcov_end(-1),
	  order_beg(0), order_end(-1),
	  plobjective_beg(0), plobjective_end(-1),
	  pldiscount_beg(0), pldiscount_end(-1)
{
	// global parse
	try {
		parse_glob(len, stream);
	} catch (const ogp::ParserException& e) {
		throw ogp::ParserException(e, dynglob_lloc.off);
	}
	// setting parameters parse
	try {
		if (paramset_end > paramset_beg)
			paramset.parse(paramset_end-paramset_beg, stream+paramset_beg);
	} catch (const ogp::ParserException& e) {
		throw ogp::ParserException(e, paramset_beg);
	}
	// model parse
	try {
		if (model_end > model_beg)
			eqs.parse(model_end-model_beg, stream+model_beg);
		else
			throw ogp::ParserException("Model section not found.", 0);
	} catch (const ogp::ParserException& e) {
		throw ogp::ParserException(e, model_beg);
	}
	// initval setting parse
	try {
		if (initval_end > initval_beg)
			initval.parse(initval_end-initval_beg, stream+initval_beg);
	} catch (const ogp::ParserException& e) {
		throw ogp::ParserException(e, initval_beg);
	}
	// vcov parse
	try {
		if (vcov_end > vcov_beg) {
			vcov.parse(vcov_end-vcov_beg, stream+vcov_beg);
		}
	} catch (const ogp::ParserException& e) {
		throw ogp::ParserException(e, vcov_beg);
	}
	// planner objective parse
	try {
		if (plobjective_end > plobjective_beg) {
			eqs.parse(plobjective_end-plobjective_beg, stream+plobjective_beg);
			t_plobjective = eqs.pop_last_formula();
		}
	} catch (const ogp::ParserException& e) {
		throw ogp::ParserException(e, plobjective_beg);
	}
	// planner discount parse
	try {
		if (pldiscount_end > pldiscount_beg) {
			t_pldiscount = parse_pldiscount(pldiscount_end - pldiscount_beg,
											stream + pldiscount_beg);
		}
	} catch (const ogp::ParserException& e) {
		throw ogp::ParserException(e, pldiscount_beg);
	}
	// order parse
	try {
		if (order_end > order_beg) {
			order = parse_order(order_end > order_beg, stream + order_beg);
		}
	} catch (const ogp::ParserException& e) {
		throw ogp::ParserException(e, order_beg);
	}

	// check the overridden order
	if (ord != -1)
		order = ord;

	// end parsing job, add planner's FOCs, make substitutions
	DynareModel::final_job();

	// calculate parameters
	calc_params();
	// calculate initial values
	calc_init();

	if (vcov_end > vcov_beg)
		vcov_mat = new ParsedMatrix(vcov);
	else {
		// vcov has not been asserted, set it to unit matrix
		vcov_mat = new TwoDMatrix(atoms.nexo(), atoms.nexo());
		vcov_mat->unit();
	}

	// check the model
	check_model();

	// differentiate
	if (order >= 1)
		eqs.differentiate(order);
}

DynareParser::DynareParser(const DynareParser& dp)
	: DynareModel(dp),
	  pa_atoms(dp.pa_atoms), paramset(dp.paramset, pa_atoms),
	  ia_atoms(dp.ia_atoms), initval(dp.initval, ia_atoms), vcov(dp.vcov),
	  model_beg(dp.model_beg), model_end(dp.model_end),
	  paramset_beg(dp.paramset_beg), paramset_end(dp.paramset_end),
	  initval_beg(dp.initval_beg), initval_end(dp.initval_end),
	  vcov_beg(dp.vcov_beg), vcov_end(dp.vcov_end),
	  order_beg(dp.order_beg), order_end(dp.order_end),
	  plobjective_beg(dp.plobjective_beg), plobjective_end(dp.plobjective_end),
	  pldiscount_beg(dp.pldiscount_beg), pldiscount_end(dp.pldiscount_end)
{
}

DynareParser::~DynareParser()
{
}

void DynareParser::add_name(const char* name, int flag)
{
	DynareModel::add_name(name, flag);
	// register with static atoms used for atom assignements
	if (flag == 1) {
		// endogenous
		ia_atoms.register_name(name);
	} else if (flag == 2) {
		// exogenous
		ia_atoms.register_name(name);
	} else if (flag == 3) {
		// parameter
		pa_atoms.register_name(name);
		ia_atoms.register_name(name);
	} else {
		throw DynareException(__FILE__, __LINE__,
							  "Unrecognized flag value.");
	}
}

void DynareParser::error(const char* mes)
{
	// throwing zero offset since this exception will be caugth at
	// constructor
	throw ogp::ParserException(mes, 0);
}

void DynareParser::print() const
{
	DynareModel::print();
	printf("parameter atoms:\n");
	paramset.print();
	printf("initval atoms:\n");
	initval.print();
	printf("model position: %d %d\n", model_beg, model_end);
	printf("paramset position: %d %d\n", paramset_beg, paramset_end);
	printf("initval position: %d %d\n", initval_beg, initval_end);
}

/** A global symbol for passing info to the DynareParser from
 * parser. */
DynareParser* dynare_parser;

/** The declarations of functions defined in dynglob_ll.cc and
 * dynglob_tab.cc generated from dynglob.lex and dynglob.y */
void* dynglob__scan_buffer(char*, size_t);
void dynglob__destroy_buffer(void*);
void dynglob_parse();
extern ogp::location_type dynglob_lloc;

void DynareParser::parse_glob(int length, const char* stream)
{
	char* buffer = new char[length+2];
	strncpy(buffer, stream, length);
	buffer[length] = '\0';
	buffer[length+1] = '\0';
	void* p = dynglob__scan_buffer(buffer, (unsigned int)length+2);
	dynare_parser = this;
	dynglob_parse();
	delete [] buffer;
	dynglob__destroy_buffer(p);
}


int DynareParser::parse_order(int len, const char* str)
{
	char* buf = new char[len+1];
	strncpy(buf, str, len);
	buf[len] = '\0';
	int res;
	sscanf(buf, "%d", &res);
	delete [] buf;
	return res;
}

int DynareParser::parse_pldiscount(int len, const char* str)
{
	char* buf = new char[len+1];
	strncpy(buf, str, len);
	buf[len] = '\0';
	if (! atoms.is_type(buf, DynareDynamicAtoms::param))
		throw ogp::ParserException(std::string("Name ") + buf + " is not a parameter", 0);

	int t = atoms.index(buf, 0);
	if (t == -1)
		t = eqs.add_nulary(buf);

	delete [] buf;
	return t;
}	

void DynareParser::calc_params()
{
	if (param_vals)
		delete param_vals;

	param_vals = new Vector(atoms.np());
	ogp::AtomAsgnEvaluator aae(paramset);
	aae.eval();
	for (int i = 0; i < atoms.np(); i++)
		(*param_vals)[i] = aae.get_value(atoms.get_params()[i]);

	for (unsigned int i = 0; i < atoms.get_params().size(); i++)
		if (! std::isfinite((*param_vals)[i]))
			printf("dynare++: warning: value for parameter %s is not finite\n",
				   atoms.get_params()[i]);
}

void DynareParser::calc_init()
{
    // update initval atoms assignings according to substitutions
	if (atom_substs)
		initval.apply_subst(atom_substs->get_old2new());

	// calculate the vector of initial values
	if (init_vals)
		delete init_vals;
	init_vals = new Vector(atoms.ny());
	ogp::AtomAsgnEvaluator aae(initval);
	// set parameters
	for (int ip = 0; ip < atoms.np(); ip++)
		aae.set_user_value(atoms.get_params()[ip], (*param_vals)[ip]);
	// set exogenous to zeros
	for (int ie = 0; ie < atoms.nexo(); ie++)
		aae.set_user_value(atoms.get_exovars()[ie], 0.0);
	// evaluate
	aae.eval();
	// set results to internally ordered vector init_vals
	for (int outer = 0; outer < atoms.ny(); outer++) {
		int i = atoms.outer2y_endo()[outer];
		(*init_vals)[i] = aae.get_value(atoms.get_endovars()[outer]);
	}

	// if the planner's FOCs have been added, then add estimate of
	// Lagrange multipliers to the vector
	if (pbuilder) {
		MultInitSS mis(*pbuilder, *param_vals, *init_vals);
	}

	// if forward substitution builder has been created, we have to
	// its substitutions and evaluate them
	if (fbuilder)
		ogdyn::DynareSteadySubstitutions dss(atoms, eqs.getTree(),
											 fbuilder->get_aux_map(), *param_vals, *init_vals);

	for (unsigned int i = 0; i < atoms.get_endovars().size(); i++)
		if (! std::isfinite((*init_vals)[i]))
			printf("dynare++: warning: initval for <%s> is not finite\n",
				   atoms.get_endovars()[atoms.y2outer_endo()[i]]);
}

// this returns false for linear functions
bool NLSelector::operator()(int t) const
{
	const ogp::Operation& op = model.getParser().getTree().operation(t);
	const DynareDynamicAtoms& atoms = model.getAtoms();
	// if the term is constant, return false
	if (model.is_constant_term(t))
		return false;
 	int nary = op.nary();
	if (nary == 0) {
		if (atoms.is_type(atoms.name(t), DynareDynamicAtoms::endovar) ||
			atoms.is_type(atoms.name(t), DynareDynamicAtoms::exovar))
			return true;
		else
			return false;
	} else if (nary == 1) {
		if (op.getCode() == ogp::UMINUS)
			return false;
		else
			return true;
	} else {
		if (op.getCode() == ogp::TIMES)
			// if at least one operand is constant, than the TIMES is linear
			if (model.is_constant_term(op.getOp1()) ||
				model.is_constant_term(op.getOp2()))
				return false;
			else
				return true;
			// both PLUS and MINUS are linear
		if (op.getCode() == ogp::PLUS ||
			op.getCode() == ogp::MINUS)
			return false;
		// POWER is linear if exponent or base is 0 or one
		if (op.getCode() == ogp::POWER &&
			(op.getOp1() == ogp::OperationTree::zero ||
			 op.getOp1() == ogp::OperationTree::one ||
			 op.getOp2() == ogp::OperationTree::zero ||
			 op.getOp2() == ogp::OperationTree::one))
			return false;
		else
			return true;
		// DIVIDE is linear if the denominator is constant, or if
		// the nominator is zero
		if (op.getCode() == ogp::DIVIDE &&
			(op.getOp1() == ogp::OperationTree::zero ||
			 model.is_constant_term(op.getOp2())))
			return false;
		else
			return true;
	}

	throw DynareException(__FILE__, __LINE__,
						  "Wrong operation in operation tree");
	return false;
}

DynareSPModel::DynareSPModel(const char** endo, int num_endo,
							 const char** exo, int num_exo,
							 const char** par, int num_par,
							 const char* equations, int len,
							 int ord)
	: DynareModel()
{
	// set the order
	order = ord;

	// add names
	for (int i = 0; i < num_endo; i++)
		add_name(endo[i], 1);
	for (int i = 0; i < num_exo; i++)
		add_name(exo[i], 2);
	for (int i = 0; i < num_par; i++)
		add_name(par[i], 3);

	// parse the equations
	eqs.parse(len, equations);

	// parsing finished
	atoms.parsing_finished(ogp::VarOrdering::bfspbfpb);

	// create what has to be created from DynareModel
	param_vals = new Vector(atoms.np());
	init_vals = new Vector(atoms.ny());
	vcov_mat = new TwoDMatrix(atoms.nexo(), atoms.nexo());

	// check the model
	check_model();

	// differentiate
	if (order >= 1)
		eqs.differentiate(order);
}

void ModelSSWriter::write_der0(FILE* fd)
{
	write_der0_preamble(fd);
	write_atom_assignment(fd);

	stop_set.clear();
	for (int fi = 0; fi < model.eqs.nformulas(); fi++)
		otree.print_operation_tree(model.eqs.formula(fi), fd, *this);

	write_der0_assignment(fd);
}

void ModelSSWriter::write_der1(FILE* fd)
{
	write_der1_preamble(fd);
	write_atom_assignment(fd);

	stop_set.clear();

	const vector<int>& variables = model.getAtoms().variables();
	const vector<int>& eam = model.getAtoms().get_endo_atoms_map();
	for (int i = 0; i < model.getParser().nformulas(); i++) {
		const ogp::FormulaDerivatives& fder = model.getParser().derivatives(i);
		for (unsigned int j = 0; j < eam.size(); j++) {
			int t = fder.derivative(ogp::FoldMultiIndex(variables.size(), 1, eam[j]));
			if (t > 0)
				otree.print_operation_tree(t, fd, *this);
		}
	}

	write_der1_assignment(fd);
}

MatlabSSWriter::MatlabSSWriter(const DynareModel& dm, const char* idd)
	: ModelSSWriter(dm), id(new char[strlen(idd)+1])
{
	strcpy(id, idd);
}


void MatlabSSWriter::write_der0_preamble(FILE* fd) const
{
	fprintf(fd,
			"%% Usage:\n"
			"%%       out = %s_f(params, y)\n"
			"%%   where\n"
			"%%       out    is a (%d,1) column vector of the residuals\n"
            "%%              of the static system\n",
			id, model.getAtoms().ny());
	write_common1_preamble(fd);
	fprintf(fd,
			"function out = %s_f(params, y)\n", id);
	write_common2_preamble(fd);
}

void MatlabSSWriter::write_der1_preamble(FILE* fd) const
{
	fprintf(fd,
			"%% Usage:\n"
			"%%       out = %s_ff(params, y)\n"
			"%%   where\n"
			"%%       out    is a (%d,%d) matrix of the first order\n"
			"%%              derivatives of the static system residuals\n"
			"%%              columns correspond to endo variables in\n"
            "%%              the ordering as declared\n",
			id, model.getAtoms().ny(), model.getAtoms().ny());
	write_common1_preamble(fd);
	fprintf(fd,
			"function out = %s_ff(params, y)\n", id);
	write_common2_preamble(fd);
}

void MatlabSSWriter::write_common1_preamble(FILE* fd) const
{
	fprintf(fd,
			"%%       params is a (%d,1) vector of parameter values\n"
			"%%              in the ordering as declared\n"
			"%%       y      is a (%d,1) vector of endogenous variables\n"
			"%%              in the ordering as declared\n"
			"%%\n"
			"%% Created by Dynare++ v. %s\n", model.getAtoms().np(),
			model.getAtoms().ny(), DYNVERSION);
	// write ordering of parameters
	fprintf(fd, "\n%% params ordering\n%% =====================\n");
	for (unsigned int ip = 0; ip < model.getAtoms().get_params().size(); ip++) {
		const char* parname = model.getAtoms().get_params()[ip];
		fprintf(fd, "%% %s\n", parname);
	}
	// write endogenous variables
	fprintf(fd, "%%\n%% y ordering\n%% =====================\n");
	for (unsigned int ie = 0; ie < model.getAtoms().get_endovars().size(); ie++) {
		const char* endoname = model.getAtoms().get_endovars()[ie];
		fprintf(fd, "%% %s\n", endoname);
	}
	fprintf(fd,"\n");
}

void MatlabSSWriter::write_common2_preamble(FILE* fd) const
{
	fprintf(fd, "if size(y) ~= [%d,1]\n\terror('Wrong size of y, must be [%d,1]');\nend\n",
			model.getAtoms().ny(), model.getAtoms().ny());
	fprintf(fd, "if size(params) ~= [%d,1]\n\terror('Wrong size of params, must be [%d,1]');\nend\n\n",
			model.getAtoms().np(), model.getAtoms().np());
}

void MatlabSSWriter::write_atom_assignment(FILE* fd) const
{
	// write OperationTree::num_constants
	fprintf(fd, "%% hardwired constants\n");
	ogp::EvalTree etree(model.getParser().getTree(), ogp::OperationTree::num_constants-1);
	for (int i = 0; i < ogp::OperationTree::num_constants; i++) {
		format_nulary(i, fd);
		double g = etree.eval(i);
		if (std::isnan(g))
			fprintf(fd, " = NaN;\n");
		else
			fprintf(fd, " = %12.8g;\n", etree.eval(i));		
	}
	// write numerical constants
	fprintf(fd, "%% numerical constants\n");
	const ogp::Constants::Tconstantmap& cmap = model.getAtoms().get_constantmap();
	for (ogp::Constants::Tconstantmap::const_iterator it = cmap.begin();
		 it != cmap.end(); ++it) {
		format_nulary((*it).first, fd);
		fprintf(fd, " = %12.8g;\n", (*it).second);
	}
	// write parameters
	fprintf(fd, "%% parameter values\n");
	for (unsigned int ip = 0; ip < model.getAtoms().get_params().size(); ip++) {
		const char* parname = model.getAtoms().get_params()[ip];
		int t = model.getAtoms().index(parname, 0);
		if (t == -1) {
			fprintf(fd, "%% %s not used in the model\n", parname);
		} else {
			format_nulary(t, fd);
			fprintf(fd, " = params(%d); %% %s\n", ip+1, parname);
		}
	}
	// write exogenous variables
	fprintf(fd, "%% exogenous variables to zeros\n");
	for (unsigned int ie = 0; ie < model.getAtoms().get_exovars().size(); ie++) {
		const char* exoname = model.getAtoms().get_exovars()[ie];
		try {
			const ogp::DynamicAtoms::Tlagmap& lmap = model.getAtoms().lagmap(exoname);
			for (ogp::DynamicAtoms::Tlagmap::const_iterator it = lmap.begin();
				 it != lmap.end(); ++it) {
				format_nulary((*it).second, fd);
				fprintf(fd, " = 0.0; %% %s\n", exoname);
			}
		} catch (const ogu::Exception& e) {
			// ignore the error of not found variable in the tree
		}
	}
	// write endogenous variables
	fprintf(fd, "%% endogenous variables to y\n");
	for (unsigned int ie = 0; ie < model.getAtoms().get_endovars().size(); ie++) {
		const char* endoname = model.getAtoms().get_endovars()[ie];
		const ogp::DynamicAtoms::Tlagmap& lmap = model.getAtoms().lagmap(endoname);
		for (ogp::DynamicAtoms::Tlagmap::const_iterator it = lmap.begin();
			 it != lmap.end(); ++it) {
			format_nulary((*it).second, fd);
			fprintf(fd, " = y(%d); %% %s\n", ie+1, endoname);
		}
	}
	fprintf(fd,"\n");
}

void MatlabSSWriter::write_der0_assignment(FILE* fd) const
{

	// initialize out variable
	fprintf(fd, "%% setting the output variable\n");
	fprintf(fd, "out = zeros(%d, 1);\n", model.getParser().nformulas());

	// fill out with the terms
	for (int i = 0; i < model.getParser().nformulas(); i++) {
		fprintf(fd, "out(%d) = ", i+1);
		format_term(model.getParser().formula(i), fd);
		fprintf(fd, ";\n");
	}
}

void MatlabSSWriter::write_der1_assignment(FILE* fd) const
{
	// initialize out variable
	fprintf(fd, "%% setting the output variable\n");
	fprintf(fd, "out = zeros(%d, %d);\n", model.getParser().nformulas(), model.getAtoms().ny());

	// fill out with the terms
	const vector<int>& variables = model.getAtoms().variables();
	const vector<int>& eam = model.getAtoms().get_endo_atoms_map();
	for (int i = 0; i < model.getParser().nformulas(); i++) {
		const ogp::FormulaDerivatives& fder = model.getParser().derivatives(i);
		for (unsigned int j = 0; j < eam.size(); j++) {
			int tvar = variables[eam[j]];
			const char* name = model.getAtoms().name(tvar);
			int yi = model.getAtoms().name2outer_endo(name);
			int t = fder.derivative(ogp::FoldMultiIndex(variables.size(), 1, eam[j]));
			if (t != ogp::OperationTree::zero) {
				fprintf(fd, "out(%d,%d) = out(%d,%d) + ", i+1, yi+1, i+1, yi+1);
				format_term(t, fd);
				fprintf(fd, "; %% %s(%d)\n", name, model.getAtoms().lead(tvar));
			}
		}
	}
}

void MatlabSSWriter::format_term(int t, FILE* fd) const
{
	fprintf(fd, "t%d", t);
}

void MatlabSSWriter::format_nulary(int t, FILE* fd) const
{
	fprintf(fd, "a%d", t);
}

void DebugOperationFormatter::format_nulary(int t, FILE* fd) const
{
	const DynareDynamicAtoms& a = model.getAtoms();

	if (t == ogp::OperationTree::zero)
		fprintf(fd, "0");
	else if (t == ogp::OperationTree::one)
		fprintf(fd, "1");
	else if (t == ogp::OperationTree::nan)
		fprintf(fd, "NaN");
	else if (t == ogp::OperationTree::two_over_pi)
		fprintf(fd, "2/sqrt(PI)");
	else if (a.is_constant(t))
		fprintf(fd, "%g", a.get_constant_value(t));
	else {
		int ll = a.lead(t);
		const char* name = a.name(t);
		if (ll == 0)
			fprintf(fd, "%s", name);
		else
			fprintf(fd, "%s(%d)", name, ll);
	}
}
