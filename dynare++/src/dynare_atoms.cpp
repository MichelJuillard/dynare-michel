// Copyright (C) 2006, Ondra Kamenik

// $Id: dynare_atoms.cpp 1765 2008-03-31 14:32:08Z kamenik $

#include "parser/cc/parser_exception.h"
#include "utils/cc/exception.h"

#include "dynare_atoms.h"

#include <string>
#include <cmath>

using namespace ogdyn;
using std::string;

void DynareStaticAtoms::register_name(const char* name)
{
	if (varnames.query(name))
		throw ogp::ParserException(string("The name ")+name+" is not unique.", 0);
	StaticAtoms::register_name(name);
}

int DynareStaticAtoms::check_variable(const char* name) const
{
	if (0 == varnames.query(name))
		throw ogp::ParserException(std::string("Unknown name <")+name+">", 0);
	Tvarmap::const_iterator it = vars.find(name);
	if (it == vars.end())
		return -1;
	else
		return (*it).second;
}

DynareDynamicAtoms::DynareDynamicAtoms(const DynareDynamicAtoms& dda)
	: SAtoms(dda)
{
	// fill atom_type
	for (Tatypemap::const_iterator it = dda.atom_type.begin();
		 it != dda.atom_type.end(); ++it)
		atom_type.insert(Tatypemap::value_type(varnames.query((*it).first), (*it).second));
}

void DynareDynamicAtoms::parse_variable(const char* in, std::string& out, int& ll) const
{
	ll = 0;
	std::string str = in;
	int left = str.find_first_of("({");
	if (left != -1) {
		out = str.substr(0, left);
		left++;
		int right = str.find_first_of(")}", left);
		if ((int)string::npos == right)
			throw ogp::ParserException(
				string("Syntax error when parsing Dynare atom <")+in+">.", 0);
		std::string tmp(str, left, right-left);
		sscanf(tmp.c_str(), "%d", &ll);
	} else {
		out = in;
	}
}

void DynareDynamicAtoms::register_uniq_endo(const char* name)
{
	FineAtoms::register_uniq_endo(name);
	atom_type.insert(Tatypemap::value_type(varnames.query(name), endovar));
}

void DynareDynamicAtoms::register_uniq_exo(const char* name)
{
	FineAtoms::register_uniq_exo(name);
	atom_type.insert(Tatypemap::value_type(varnames.query(name), exovar));
}

void DynareDynamicAtoms::register_uniq_param(const char* name)
{
	FineAtoms::register_uniq_param(name);
	atom_type.insert(Tatypemap::value_type(varnames.query(name), param));
}

bool DynareDynamicAtoms::is_type(const char* name, atype tp) const
{
	Tatypemap::const_iterator it = atom_type.find(name);
	if (it != atom_type.end() && (*it).second == tp)
		return true;
	else
		return false;
}

void DynareDynamicAtoms::print() const
{
	SAtoms::print();
	printf("Name types:\n");
	for (Tatypemap::const_iterator it = atom_type.begin();
		 it != atom_type.end(); ++it)
		printf("name=%s type=%s\n", (*it).first,
			   ((*it).second == endovar) ? "endovar" : (((*it).second == exovar)? "exovar" : "param"));
}

std::string DynareDynamicAtoms::convert(int t) const
{
	if (t < ogp::OperationTree::num_constants) {
		throw ogu::Exception(__FILE__,__LINE__,
							 "Tree index is a built-in constant in DynareDynamicAtoms::convert");
		return std::string();
	}
	if (is_constant(t)) {
		double v = get_constant_value(t);
		char buf[100];
		sprintf(buf, "%20.16g", v);
		const char* s = buf;
		while (*s == ' ')
			++s;
		return std::string(s);
	}
	
	const char* s = name(t);
	if (is_type(s, endovar)) {
		int ll = lead(t);
		char buf[100];
		if (ll)
			sprintf(buf, "%s(%d)", s, ll);
		else
			sprintf(buf, "%s", s);
		return std::string(buf);
	}

	return std::string(s);
}


void DynareAtomValues::setValues(ogp::EvalTree& et) const
{
	// set constants
	atoms.setValues(et);

	// set parameteres
	for (unsigned int i = 0; i < atoms.get_params().size(); i++) {
		try {
			const ogp::DynamicAtoms::Tlagmap& lmap = atoms.lagmap(atoms.get_params()[i]);
			for (ogp::DynamicAtoms::Tlagmap::const_iterator it = lmap.begin();
				 it != lmap.end(); ++it) {
				int t = (*it).second;
				et.set_nulary(t, paramvals[i]);
			}
		} catch (const ogu::Exception& e) {
			// ignore non-referenced parameters; there is no
			// lagmap for them
		}
	}

	// set endogenous
	for (unsigned int outer_i = 0; outer_i < atoms.get_endovars().size(); outer_i++) {
		try {
			const ogp::DynamicAtoms::Tlagmap& lmap = atoms.lagmap(atoms.get_endovars()[outer_i]);
			for (ogp::DynamicAtoms::Tlagmap::const_iterator it = lmap.begin();
				 it != lmap.end(); ++it) {
				int ll = (*it).first;
				int t = (*it).second;
				int i = atoms.outer2y_endo()[outer_i];
				if (ll == -1) {
					et.set_nulary(t, yym[i-atoms.nstat()]);
				}
				else if (ll == 0)
					et.set_nulary(t, yy[i]);
				else
					et.set_nulary(t, yyp[i-atoms.nstat()-atoms.npred()]);
			}
		} catch (const ogu::Exception& e) {
			// ignore non-referenced endogenous variables; there is no
			// lagmap for them
		}
	}

	// set exogenous
	for (unsigned int outer_i = 0; outer_i < atoms.get_exovars().size(); outer_i++) {
		try {
			const ogp::DynamicAtoms::Tlagmap& lmap = atoms.lagmap(atoms.get_exovars()[outer_i]);
			for (ogp::DynamicAtoms::Tlagmap::const_iterator it = lmap.begin();
				 it != lmap.end(); ++it) {
				int ll = (*it).first;
				if (ll == 0) { // this is always true because of checks
					int t = (*it).second;
					int i = atoms.outer2y_exo()[outer_i];			
					et.set_nulary(t, xx[i]);
				}
			}
		} catch (const ogu::Exception& e) {
			// ignore non-referenced variables
		}
	}
}

void DynareStaticSteadyAtomValues::setValues(ogp::EvalTree& et) const
{
	// set constants
	atoms_static.setValues(et);

	// set parameters
	for (unsigned int i = 0; i < atoms_static.get_params().size(); i++) {
		const char* name = atoms_static.get_params()[i];
		int t = atoms_static.index(name);
		if (t != -1) {
			int idyn = atoms.name2outer_param(name);
			et.set_nulary(t, paramvals[idyn]);
		}
	}

	// set endogenous
	for (unsigned int i = 0; i < atoms_static.get_endovars().size(); i++) {
		const char* name = atoms_static.get_endovars()[i];
		int t = atoms_static.index(name);
		if (t != -1) {
			int idyn = atoms.outer2y_endo()[atoms.name2outer_endo(name)];
			et.set_nulary(t, yy[idyn]);
		}
	}

	// set exogenous
	for (unsigned int i = 0; i < atoms_static.get_exovars().size(); i++) {
		const char* name = atoms_static.get_exovars()[i];
		int t = atoms_static.index(name);
		if (t != -1)
			et.set_nulary(t, 0.0);
	}
}

DynareSteadySubstitutions::DynareSteadySubstitutions(const ogp::FineAtoms& a,
													 const ogp::OperationTree& tree,
													 const Tsubstmap& subst,
													 const Vector& pvals, Vector& yy)
	: atoms(a), y(yy)
{
	// fill the vector of left and right hand sides
	for (Tsubstmap::const_iterator it = subst.begin();
		 it != subst.end(); ++it) {
		left_hand_sides.push_back((*it).first);
		right_hand_sides.push_back((*it).second);
	}

	// evaluate right hand sides
	DynareSteadyAtomValues dsav(atoms, pvals, y);
	ogp::FormulaCustomEvaluator fe(tree, right_hand_sides);
	fe.eval(dsav, *this);
}

void DynareSteadySubstitutions::load(int i, double res)
{
	const char* name = left_hand_sides[i];
	int iouter = atoms.name2outer_endo(name);
	int iy = atoms.outer2y_endo()[iouter];
	if (! std::isfinite(y[iy]))
		y[iy] = res;
}

DynareStaticSteadySubstitutions::
DynareStaticSteadySubstitutions(const ogp::FineAtoms& a, const ogp::StaticFineAtoms& sa,
								const ogp::OperationTree& tree,
								const Tsubstmap& subst,
								const Vector& pvals, Vector& yy)
	: atoms(a), atoms_static(sa), y(yy)
{
	// fill the vector of left and right hand sides
	for (Tsubstmap::const_iterator it = subst.begin();
		 it != subst.end(); ++it) {
		left_hand_sides.push_back((*it).first);
		right_hand_sides.push_back((*it).second);
	}

	// evaluate right hand sides
	DynareStaticSteadyAtomValues dsav(atoms, atoms_static, pvals, y);
	ogp::FormulaCustomEvaluator fe(tree, right_hand_sides);
	fe.eval(dsav, *this);
}

void DynareStaticSteadySubstitutions::load(int i, double res)
{
	const char* name = left_hand_sides[i];
	int iouter = atoms.name2outer_endo(name);
	int iy = atoms.outer2y_endo()[iouter];
	if (! std::isfinite(y[iy]))
		y[iy] = res;
}
