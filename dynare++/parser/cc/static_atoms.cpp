// Copyright (C) 2006, Ondra Kamenik

// $Id: static_atoms.cpp 1360 2007-07-10 11:44:20Z kamenik $

#include "static_atoms.h"
#include "utils/cc/exception.h"

using namespace ogp;

StaticAtoms::StaticAtoms(const StaticAtoms& a)
	: Atoms(), Constants(a), varnames(a.varnames),
	  varorder(), vars(), indices()
{
	// fill varorder
	for (unsigned int i = 0; i < a.varorder.size(); i++) {
		const char* s = varnames.query(a.varorder[i]);
		varorder.push_back(s);
	}

	// fill vars
	for (Tvarmap::const_iterator it = a.vars.begin();
		 it != a.vars.end(); ++it) {
		const char* s = varnames.query((*it).first);
		vars.insert(Tvarmap::value_type(s, (*it).second));
	}

	// fill indices
	for (Tinvmap::const_iterator it = a.indices.begin();
		 it != a.indices.end(); ++it) {
		const char* s = varnames.query((*it).second);
		indices.insert(Tinvmap::value_type((*it).first, s));
	}
}

void StaticAtoms::import_atoms(const DynamicAtoms& da, OperationTree& otree, Tintintmap& tmap)
{
	Constants::import_constants(da, otree, tmap);

	for (int i = 0; i < da.get_name_storage().num(); i++) {
		const char* name = da.get_name_storage().get_name(i);
		register_name(name);
		int tnew = otree.add_nulary();
		assign(name, tnew);
		try {
			const DynamicAtoms::Tlagmap& lmap = da.lagmap(name);
			for (DynamicAtoms::Tlagmap::const_iterator it = lmap.begin();
				 it != lmap.end(); ++it) {
				int told = (*it).second;
				tmap.insert(Tintintmap::value_type(told, tnew));
			}
		} catch (const ogu::Exception& e) {
		}
	}
}


int StaticAtoms::check(const char* name) const
{
	if (DynamicAtoms::is_string_constant(name)) {
		return Constants::check(name);
	} else {
		return check_variable(name);
	}
}

int StaticAtoms::index(const char* name) const
{
	Tvarmap::const_iterator it = vars.find(name);
	if (it == vars.end())
		return -1;
	else
		return (*it).second;
}

const char* StaticAtoms::inv_index(int t) const
{
	Tinvmap::const_iterator it = indices.find(t);
	if (it == indices.end())
		return NULL;
	else
		return (*it).second;
}

void StaticAtoms::assign(const char* name, int t)
{
	if (DynamicAtoms::is_string_constant(name)) {
		double val;
		sscanf(name, "%lf", &val);
		add_constant(t, val);
	} else {
		const char* ss = varnames.insert(name);
		vars.insert(Tvarmap::value_type(ss, t));
		indices.insert(Tinvmap::value_type(t, ss));
	}
}

vector<int> StaticAtoms::variables() const
{
	vector<int> res;
	for (Tvarmap::const_iterator it = vars.begin();
		 it != vars.end(); ++it) {
		res.push_back((*it).second);
	}
	return res;
}

void StaticAtoms::register_name(const char* name)
{
	const char* ss = varnames.insert(name);
	varorder.push_back(ss);
}

void StaticAtoms::print() const
{
	printf("constants:\n");
	Constants::print();
	printf("variable names:\n");
	varnames.print();
	printf("map to tree indices:\n");
	for (Tvarmap::const_iterator it = vars.begin(); it != vars.end(); ++it)
		printf("%s\t->\t%d\n", (*it).first, (*it).second);
}
