// Copyright (C) 2006, Ondra Kamenik

// $Id: static_fine_atoms.cpp 82 2007-04-19 11:33:30Z ondra $

#include "utils/cc/exception.h"

#include "static_fine_atoms.h"
#include "parser_exception.h"

using namespace ogp;

StaticFineAtoms::StaticFineAtoms(const StaticFineAtoms& sfa)
	: StaticAtoms(sfa),
	  params(), param_outer_map(),
	  endovars(), endo_outer_map(),
	  exovars(), exo_outer_map(),
	  der_atoms(sfa.der_atoms),
	  endo_atoms_map(sfa.endo_atoms_map),
	  exo_atoms_map(sfa.exo_atoms_map)
{
	for (unsigned int i = 0; i < sfa.params.size(); i++) {
		const char* name = varnames.query(sfa.params[i]);
		params.push_back(name);
		param_outer_map.insert(Tvarintmap::value_type(name, i));
	}

	for (unsigned int i = 0; i < sfa.endovars.size(); i++) {
		const char* name = varnames.query(sfa.endovars[i]);
		endovars.push_back(name);
		endo_outer_map.insert(Tvarintmap::value_type(name, i));
	}

	for (unsigned int i = 0; i < sfa.exovars.size(); i++) {
		const char* name = varnames.query(sfa.exovars[i]);
		exovars.push_back(name);
		exo_outer_map.insert(Tvarintmap::value_type(name, i));
	}
}

void StaticFineAtoms::import_atoms(const FineAtoms& fa, OperationTree& otree, Tintintmap& tmap)
{
	StaticAtoms::import_atoms(fa, otree, tmap);

	// we just need to put parameters, endovars, and exovars to
	// respective vectors, the names are already in the storage

	// parameters
	const vector<const char*>& fa_params = fa.get_params();
	for (unsigned int i = 0; i < fa_params.size(); i++)
		register_param(fa_params[i]);

	// endogenous
	const vector<const char*>& fa_endovars = fa.get_endovars();
	for (unsigned int i = 0; i < fa_endovars.size(); i++)
		register_endo(fa_endovars[i]);

	// exogenous
	const vector<const char*>& fa_exovars = fa.get_exovars();
	for (unsigned int i = 0; i < fa_exovars.size(); i++)
		register_exo(fa_exovars[i]);

	parsing_finished();
}

void StaticFineAtoms::import_atoms(const FineAtoms& fa, OperationTree& otree, Tintintmap& tmap,
								   const char* dummy)
{
	StaticAtoms::import_atoms(fa, otree, tmap);

	// we just need to put parameters, endovars, and exovars to
	// respective vectors, the names are already in the storage

	// parameters
	const vector<const char*>& fa_params = fa.get_params();
	for (unsigned int i = 0; i < fa_params.size(); i++)
		register_param(fa_params[i]);

	// endogenous
	const vector<const char*>& fa_endovars = fa.get_endovars();
	for (unsigned int i = 0; i < fa_endovars.size(); i++)
		register_endo(fa_endovars[fa.y2outer_endo()[i]]);

	// exogenous
	const vector<const char*>& fa_exovars = fa.get_exovars();
	for (unsigned int i = 0; i < fa_exovars.size(); i++)
		register_exo(fa_exovars[fa.y2outer_exo()[i]]);

	parsing_finished();
}

int StaticFineAtoms::check_variable(const char* name) const
{
	const char* ss = varnames.query(name);
	if (ss == NULL)
		throw ParserException(string("Variable <")+name+"> not declared.",0);
	return index(name);
}

void StaticFineAtoms::parsing_finished()
{
	// build der_atoms, and endo_atoms_map and exo_atoms_map
	der_atoms.clear();
	endo_atoms_map.clear();
	exo_atoms_map.clear();

	// go through all endo and exo insert tree indices, ignore names
	// whose tree index is -1 (those which are not referenced)
	for (unsigned int i = 0; i < endovars.size(); i++) {
		int t = index(endovars[i]);
		if (t != -1) {
			endo_atoms_map.push_back(der_atoms.size());
			der_atoms.push_back(t);
		}
	}
	for (unsigned int i = 0; i < exovars.size(); i++) {
		int t = index(exovars[i]);
		if (t != -1) {
			exo_atoms_map.push_back(der_atoms.size());
			der_atoms.push_back(t);
		}
	}
}

int StaticFineAtoms::name2outer_param(const char* name) const
{
	Tvarintmap::const_iterator it = param_outer_map.find(name);
	if (it == param_outer_map.end())
		throw ogu::Exception(__FILE__,__LINE__,
							 "Name is not a parameter in StaticFineAtoms::name2outer_param");
	return (*it).second;
}

int StaticFineAtoms::name2outer_endo(const char* name) const
{
	Tvarintmap::const_iterator it = endo_outer_map.find(name);
	if (it == endo_outer_map.end())
		throw ogu::Exception(__FILE__,__LINE__,
							 "Name is not an endogenous variable in StaticFineAtoms::name2outer_endo");
	return (*it).second;
}

int StaticFineAtoms::name2outer_exo(const char* name) const
{
	Tvarintmap::const_iterator it = exo_outer_map.find(name);
	if (it == exo_outer_map.end())
		throw ogu::Exception(__FILE__,__LINE__,
							 "Name is not an exogenous variable in StaticFineAtoms::name2outer_exo");
	return (*it).second;
}

void StaticFineAtoms::register_uniq_endo(const char* name)
{
	if (varnames.query(name))
		throw ogp::ParserException(string("Endogenous variable <")+name+"> is not unique.",0);
	const char* ss = varnames.insert(name);
	register_endo(ss);
}

void StaticFineAtoms::register_uniq_exo(const char* name)
{
	if (varnames.query(name))
		throw ogp::ParserException(string("Exogenous variable <")+name+"> is not unique.",0);
	const char* ss = varnames.insert(name);
	register_exo(ss);
}

void StaticFineAtoms::register_uniq_param(const char* name)
{
	if (varnames.query(name))
		throw ogp::ParserException(string("Parameter <")+name+"> is not unique.",0);
	const char* ss = varnames.insert(name);
	register_param(ss);
}

void StaticFineAtoms::print() const
{
	StaticAtoms::print();
	printf("endo atoms map:\n");
	for (unsigned int i = 0; i < endo_atoms_map.size(); i++)
		printf("%d --> %d\n", i, endo_atoms_map[i]);
	printf("exo atoms map:\n");
	for (unsigned int i = 0; i < exo_atoms_map.size(); i++)
		printf("%d --> %d\n", i, exo_atoms_map[i]);	
	printf("der atoms:\n");
	for (unsigned int i = 0; i < der_atoms.size(); i++)
		printf("%d\t%d\n",i, der_atoms[i]);
}

void StaticFineAtoms::register_endo(const char* name)
{
	const char* ss = varnames.query(name);
	if (ss == NULL)
		throw ogp::ParserException(string("Endogenous variable <")
								   +name+"> not found in storage.",0);
	endovars.push_back(ss);
	endo_outer_map.insert(Tvarintmap::value_type(ss, endovars.size()-1));
}

void StaticFineAtoms::register_exo(const char* name)
{
	const char* ss = varnames.query(name);
	if (ss == NULL)
		throw ogp::ParserException(string("Exogenous variable <")
								   +name+"> not found in storage.",0);
	exovars.push_back(ss);
	exo_outer_map.insert(Tvarintmap::value_type(ss, exovars.size()-1));
}

void StaticFineAtoms::register_param(const char* name)
{
	const char* ss = varnames.query(name);
	if (ss == NULL)
		throw ogp::ParserException(string("Parameter <")+name+"> not found in storage.",0);
	params.push_back(ss);
	param_outer_map.insert(Tvarintmap::value_type(ss, params.size()-1));
}

