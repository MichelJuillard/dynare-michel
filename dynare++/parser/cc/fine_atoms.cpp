// Copyright (C) 2005, Ondra Kamenik

// $Id: fine_atoms.cpp 1759 2008-03-31 14:25:20Z kamenik $

#include "utils/cc/exception.h"

#include "parser_exception.h"
#include "fine_atoms.h"

using namespace ogp;

AllvarOuterOrdering::AllvarOuterOrdering(const vector<const char*>& allvar_outer,
										 const FineAtoms& a)
	: atoms(a), allvar(),
	  endo2all(a.get_endovars().size(), -1),
	  exo2all(a.get_exovars().size(), -1)
{
	// fill in the allvar from allvar_outer
	for (unsigned int i = 0; i < allvar_outer.size(); i++) {
		const char* s = atoms.varnames.query(allvar_outer[i]);
		if (s)
			allvar.push_back(s);
		else
			throw ogu::Exception(__FILE__, __LINE__,
								 string("Variable ") + allvar_outer[i] + " is not a declared symbol in AllvarOuterOrdering constructor");
	}

	// fill in endo2all and exo2all
	for (unsigned int i = 0; i < allvar.size(); i++) {
		Tvarintmap::const_iterator it = atoms.endo_outer_map.find(allvar[i]);
		if (it != atoms.endo_outer_map.end())
			endo2all[(*it).second] = i;
		else {
			it = atoms.exo_outer_map.find(allvar[i]);
			if (it != atoms.exo_outer_map.end())
				exo2all[(*it).second] = i;
			else
				throw ogu::Exception(__FILE__, __LINE__,
									 string("Name ") + allvar[i] + " is neither endogenous nor exogenous variable in AllvarOuterOrdering constructor");
		}
	}

	// check whether everything has been filled
	unsigned int iendo = 0;
	while (iendo < endo2all.size() && endo2all[iendo] != -1) iendo++;
	unsigned int iexo = 0;
	while (iexo < exo2all.size() && exo2all[iexo] != -1) iexo++;
	if (iendo < endo2all.size())
		throw ogu::Exception(__FILE__, __LINE__,
							 string("Endogenous variable ") + atoms.get_endovars()[iendo] +
							 " not found in outer all ordering in AllvarOuterOrdering constructor");
	if (iexo < exo2all.size())
		throw ogu::Exception(__FILE__, __LINE__,
							 string("Exogenous variable ") + atoms.get_exovars()[iexo] +
							 " not found in outer all ordering in AllvarOuterOrdering constructor");
}

AllvarOuterOrdering::AllvarOuterOrdering(const AllvarOuterOrdering& avo,
										 const FineAtoms& a)
	: atoms(a), allvar(),
	  endo2all(avo.endo2all),
	  exo2all(avo.exo2all)
{
	// fill in the allvar from avo.allvar
	for (unsigned int i = 0; i < avo.allvar.size(); i++) {
		const char* s = atoms.varnames.query(avo.allvar[i]);
		allvar.push_back(s);
	}
}


FineAtoms::FineAtoms(const FineAtoms& fa)
	: DynamicAtoms(fa), params(), endovars(), exovars(),
	  endo_order(NULL), exo_order(NULL), allvar_order(NULL),
	  der_atoms(fa.der_atoms),
	  endo_atoms_map(fa.endo_atoms_map),
	  exo_atoms_map(fa.exo_atoms_map)
{
	// fill in params
	for (unsigned int i = 0; i < fa.params.size(); i++) {
		const char* s = varnames.query(fa.params[i]);
		if (! s)
			throw ogu::Exception(__FILE__, __LINE__,
								 string("Parameter ") + fa.params[i] + " does not exist in FineAtoms copy cosntructor");
		params.push_back(s);
		param_outer_map.insert(Tvarintmap::value_type(s, params.size()-1));
	}
	// fill in endovars
	for (unsigned int i = 0; i < fa.endovars.size(); i++) {
		const char* s = varnames.query(fa.endovars[i]);
		if (! s)
			throw ogu::Exception(__FILE__, __LINE__,
								 string("Endo variable ") + fa.endovars[i] + " does not exist in FineAtoms copy constructor");
		endovars.push_back(s);
		endo_outer_map.insert(Tvarintmap::value_type(s, endovars.size()-1));
	}
	// fill in exovars
	for (unsigned int i = 0; i < fa.exovars.size(); i++) {
		const char* s = varnames.query(fa.exovars[i]);
		if (! s)
			throw ogu::Exception(__FILE__, __LINE__,
								 string("Exo variable ") + fa.exovars[i] + " does not exist in FineAtoms copy cosntructor");
		exovars.push_back(s);
		exo_outer_map.insert(Tvarintmap::value_type(s, exovars.size()-1));
	}

	if (fa.endo_order)
		endo_order = fa.endo_order->clone(endovars, *this);

	if (fa.exo_order)
		exo_order = fa.exo_order->clone(exovars, *this);

	if (fa.allvar_order)
		allvar_order = new AllvarOuterOrdering(*(fa.allvar_order), *this);
}

int FineAtoms::check_variable(const char* name) const
{
	string str;
	int ll;
	parse_variable(name, str, ll);
	if (varnames.query(str.c_str()))
		return DynamicAtoms::check_variable(name);
	else {
		throw ParserException(string("Variable <")+str+"> not declared.",0);
		return -1;
	}
}

int FineAtoms::num_exo_periods() const
{
	int mlead, mlag;
	exovarspan(mlead, mlag);
	return mlead-mlag+1;
}

void FineAtoms::parsing_finished(VarOrdering::ord_type ot)
{
	make_internal_orderings(ot);

	// by default, concatenate outer endo and outer exo and make it as
	// allvar outer:
	vector<const char*> allvar_tmp;
	allvar_tmp.insert(allvar_tmp.end(), endovars.begin(), endovars.end());
	allvar_tmp.insert(allvar_tmp.end(), exovars.begin(), exovars.end());

	if (allvar_order)
		delete allvar_order;
	allvar_order = new AllvarOuterOrdering(allvar_tmp, *this);
}

void FineAtoms::parsing_finished(VarOrdering::ord_type ot,
								 const vector<const char*> allvar)
{
	make_internal_orderings(ot);
	if (allvar_order)
		delete allvar_order;
	allvar_order = new AllvarOuterOrdering(allvar, *this);
}

const vector<const char*>& FineAtoms::get_allvar() const
{
	if (! allvar_order)
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::get_allvars called before parsing_finished");

	return allvar_order->get_allvar();
}

const vector<int>& FineAtoms::outer_endo2all() const
{
	if (! allvar_order)
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::outer_endo2all called before parsing_finished");

	return allvar_order->get_endo2all();
}

const vector<int>& FineAtoms::outer_exo2all() const
{
	if (! allvar_order)
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::outer_exo2all called before parsing_finished");

	return allvar_order->get_exo2all();
}


vector<int> FineAtoms::variables() const
{
	if (endo_order) {
		return der_atoms;
	} else {
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::variables called before parsing_finished");
		return vector<int>();
	}
}

int FineAtoms::nstat() const
{
	if (endo_order) {
		return endo_order->nstat();
	} else {
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::nstat called before parsing_finished");
		return -1;
	}
}

int FineAtoms::npred() const
{
	if (endo_order) {
		return endo_order->npred();
	} else {
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::npred called before parsing_finished");
		return -1;
	}
}

int FineAtoms::nboth() const
{
	if (endo_order) {
		return endo_order->nboth();
	} else {
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::nboth called before parsing_finished");
		return -1;
	}
}

int FineAtoms::nforw() const
{
	if (endo_order) {
		return endo_order->nforw();
	} else {
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::nforw called before parsing_finished");
		return -1;
	}
}

int FineAtoms::get_pos_of_endo(int t) const
{
	if (endo_order) {
		return endo_order->get_pos_of(t);
	} else {
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::get_pos_of_endo called before parsing_finished");
		return -1;
	}
}

int FineAtoms::get_pos_of_exo(int t) const
{
	if (exo_order) {
		return exo_order->get_pos_of(t);
	} else {
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::get_pos_of_exo called before parsing_finished");
		return -1;
	}
}

int FineAtoms::get_pos_of_all(int t) const
{
	if (endo_order && exo_order) {
		if (endo_order->check(t))
			return endo_order->get_pos_of(t);
		else if (exo_order->check(t))
			return endo_order->length() + exo_order->get_pos_of(t);
		else {
			throw ogu::Exception(__FILE__,__LINE__,
								 "Atom is not endo nor exo in FineAtoms::get_pos_of_all");
			return -1;
		}
	} else {
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::get_pos_of_exo called before parsing_finished");
		return -1;
	}
}

const vector<int>& FineAtoms::y2outer_endo() const
{
	if (! endo_order)
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::y2outer_endo called before parsing_finished");
	return endo_order->get_y2outer();
}

const vector<int>& FineAtoms::outer2y_endo() const
{
	if (! endo_order)
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::outer2y_endo called before parsing_finished");
	return endo_order->get_outer2y();
}

const vector<int>& FineAtoms::y2outer_exo() const
{
	if (! exo_order)
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::y2outer_endo called before parsing_finished");
	return exo_order->get_y2outer();
}

const vector<int>& FineAtoms::outer2y_exo() const
{
	if (! exo_order)
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::outer2y_exo called before parsing_finished");
	return exo_order->get_outer2y();
}

const vector<int>& FineAtoms::get_endo_atoms_map() const
{
	if (! endo_order)
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::get_endo_atoms_map called before parsing_finished");
	return endo_atoms_map;
}

const vector<int>& FineAtoms::get_exo_atoms_map() const
{
	if (! exo_order)
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::get_exo_atoms_map called before parsing_finished");
	return exo_atoms_map;
}

int FineAtoms::name2outer_param(const char* name) const
{
	Tvarintmap::const_iterator it = param_outer_map.find(name);
	if (it == param_outer_map.end())
		throw ogu::Exception(__FILE__,__LINE__,
							 "Name is not a parameter in FineAtoms::name2outer_param");
	return (*it).second;
}

int FineAtoms::name2outer_endo(const char* name) const
{
	Tvarintmap::const_iterator it = endo_outer_map.find(name);
	if (it == endo_outer_map.end())
		throw ogu::Exception(__FILE__,__LINE__,
							 "Name is not an endogenous variable in FineAtoms::name2outer_endo");
	return (*it).second;
}

int FineAtoms::name2outer_exo(const char* name) const
{
	Tvarintmap::const_iterator it = exo_outer_map.find(name);
	if (it == exo_outer_map.end())
		throw ogu::Exception(__FILE__,__LINE__,
							 "Name is not an exogenous variable in FineAtoms::name2outer_exo");
	return (*it).second;
}

int FineAtoms::name2outer_allvar(const char* name) const
{
	if (! allvar_order)
		throw ogu::Exception(__FILE__,__LINE__,
							 "FineAtoms::name2outer_allvar called beore parsing_finished");

	Tvarintmap::const_iterator it = endo_outer_map.find(name);
	if (it != endo_outer_map.end())
		return allvar_order->get_endo2all()[(*it).second];
	else {
		it = exo_outer_map.find(name);
		if (it != exo_outer_map.end())
			return allvar_order->get_exo2all()[(*it).second];
	}

	throw ogu::Exception(__FILE__,__LINE__,
						 string("Name ") + name + " is neither endo nor exo variable in FineAtoms::name2outer_allvar");
	return -1;
}

void FineAtoms::register_uniq_endo(const char* name)
{
	if (varnames.query(name))
		throw ogp::ParserException(string("Endogenous variable <")+name+"> is not unique.",0);
	const char* ss = varnames.insert(name);
	endovars.push_back(ss);
	endo_outer_map.insert(Tvarintmap::value_type(ss, endovars.size()-1));
}

void FineAtoms::register_uniq_exo(const char* name)
{
	if (varnames.query(name))
		throw ogp::ParserException(string("Exogenous variable <")+name+"> is not unique.",0);
	const char* ss = varnames.insert(name);
	exovars.push_back(ss);
	exo_outer_map.insert(Tvarintmap::value_type(ss, exovars.size()-1));
}

void FineAtoms::register_uniq_param(const char* name)
{
	if (varnames.query(name))
		throw ogp::ParserException(string("Parameter <")+name+"> is not unique.",0);
	const char* ss = varnames.insert(name);
	params.push_back(ss);
	param_outer_map.insert(Tvarintmap::value_type(ss, params.size()-1));
}

void FineAtoms::make_internal_orderings(VarOrdering::ord_type ot)
{
	bool endo_ordering_done = false;
	bool exo_ordering_done = false;

	order_type = ot;

	int mlead, mlag;
	endovarspan(mlead, mlag);
	if (mlag >= -1 && mlead <= 1) {
		// make endo ordering
		if (endo_order)
			delete endo_order;
		if (ot == VarOrdering::pbspbfbf)
			endo_order = new EndoVarOrdering1(endovars, *this);
		else
			endo_order = new EndoVarOrdering2(endovars, *this);
		endo_order->do_ordering();
		endo_ordering_done = true;
	}

	exovarspan(mlead, mlag);
	if (mlag == 0 && mlead == 0) {
		// make exo ordering
		if (exo_order)
			delete exo_order;
		exo_order = new ExoVarOrdering(exovars, *this);
		exo_order->do_ordering();
		exo_ordering_done = true;
	}

	if (endo_ordering_done && exo_ordering_done) {
		// concatenate der atoms from endo_order and exo_order
		der_atoms.clear();
		der_atoms.insert(der_atoms.end(), 
						 endo_order->get_der_atoms().begin(),
						 endo_order->get_der_atoms().end());
		der_atoms.insert(der_atoms.end(), 
						 exo_order->get_der_atoms().begin(),
						 exo_order->get_der_atoms().end());
		
		// create endo_atoms_map; der_atoms is a concatenation, so it is easy
		int endo_atoms = endo_order->get_der_atoms().size();
		endo_atoms_map.clear();
		for (int i = 0; i < endo_atoms; i++)
			endo_atoms_map.push_back(i);
		// create exo_atoms_map
		int exo_atoms = exo_order->get_der_atoms().size();
		exo_atoms_map.clear();
		for (int i = 0; i < exo_atoms; i++)
			exo_atoms_map.push_back(endo_atoms + i);
	}
}

void FineAtoms::print() const
{
	DynamicAtoms::print();
	if (endo_order) {
		printf("Endo ordering:\n");
		endo_order->print();
	} else {
		printf("Endo ordering not created.\n");
	}
	if (exo_order) {
		printf("Exo ordering:\n");
		exo_order->print();
	} else {
		printf("Exo ordering not created.\n");
	}
	printf("endo atoms map:\n");
	for (unsigned int i = 0; i < endo_atoms_map.size(); i++)
		printf("%d --> %d\n", i, endo_atoms_map[i]);
	printf("exo atoms map:\n");
	for (unsigned int i = 0; i < exo_atoms_map.size(); i++)
		printf("%d --> %d\n", i, exo_atoms_map[i]);
}
