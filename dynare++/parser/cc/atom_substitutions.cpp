// Copyright (C) 2006, Ondra Kamenik

// $Id: atom_substitutions.cpp 42 2007-01-22 21:53:24Z ondra $

#include "atom_substitutions.h"
#include "utils/cc/exception.h"

using namespace ogp;

AtomSubstitutions::AtomSubstitutions(const AtomSubstitutions& as, const FineAtoms& oa,
									 FineAtoms& na)
	: old_atoms(oa), new_atoms(na)
{
	const NameStorage& ns = na.get_name_storage();

	// fill new2old
	for (Tshiftmap::const_iterator it = as.new2old.begin();
		 it != as.new2old.end(); ++it)
		new2old.insert(Tshiftmap::value_type(ns.query((*it).first),
											 Tshiftname(ns.query((*it).second.first),
														(*it).second.second)));
	// fill old2new
	for (Toldnamemap::const_iterator it = as.old2new.begin();
		 it != as.old2new.end(); ++it) {
		Tshiftnameset sset;
		for (Tshiftnameset::const_iterator itt = (*it).second.begin();
			 itt != (*it).second.end(); ++itt)
			sset.insert(Tshiftname(ns.query((*itt).first), (*itt).second));
		old2new.insert(Toldnamemap::value_type(ns.query((*it).first), sset));
	}
}


void AtomSubstitutions::add_substitution(const char* newname, const char* oldname, int tshift)
{
	// make sure the storage is from the new_atoms
	newname = new_atoms.get_name_storage().query(newname);
	oldname = new_atoms.get_name_storage().query(oldname);
	if (newname == NULL || oldname == NULL)
		throw ogu::Exception(__FILE__,__LINE__,
							 "Bad newname or oldname in AtomSubstitutions::add_substitution");

	// insert to new2old map
	new2old.insert(Tshiftmap::value_type(newname, Tshiftname(oldname, tshift)));
	// insert to old2new map
	Toldnamemap::iterator it = old2new.find(oldname);
	if (it != old2new.end())
		(*it).second.insert(Tshiftname(newname, -tshift));
	else {
		Tshiftnameset snset;
		snset.insert(Tshiftname(newname, -tshift));
		old2new.insert(Toldnamemap::value_type(oldname, snset));
	}

	// put to info
	info.num_substs++;
}

void AtomSubstitutions::substitutions_finished(VarOrdering::ord_type ot)
{
	// create an external ordering of new_atoms from old_atoms
	const vector<const char*>& oa_ext = old_atoms.get_allvar();
	vector<const char*> na_ext;
	for (unsigned int i = 0; i < oa_ext.size(); i++) {
		const char* oname = oa_ext[i];
		// add the old name itself
		na_ext.push_back(oname);
		// add all new names derived from the old name
		Toldnamemap::const_iterator it = old2new.find(oname);
		if (it != old2new.end())
			for (Tshiftnameset::const_iterator itt = (*it).second.begin();
				 itt != (*it).second.end(); ++itt)
				na_ext.push_back((*itt).first);
	}

	// call parsing finished for the new_atoms
	new_atoms.parsing_finished(ot, na_ext);
}

const char* AtomSubstitutions::get_new4old(const char* oldname, int tshift) const
{
	Toldnamemap::const_iterator it = old2new.find(oldname);
	if (it != old2new.end()) {
		const Tshiftnameset& sset = (*it).second;
		for (Tshiftnameset::const_iterator itt = sset.begin();
			 itt != sset.end(); ++itt)
			if ((*itt).second == - tshift)
				return (*itt).first;
	}
	return NULL;
}

void AtomSubstitutions::print() const
{
	printf("Atom Substitutions:\nOld ==> New:\n");
	for (Toldnamemap::const_iterator it = old2new.begin(); it != old2new.end(); ++it)
		for (Tshiftnameset::const_iterator itt = (*it).second.begin();
			 itt != (*it).second.end(); ++itt)
			printf("    %s ==> [%s, %d]\n", (*it).first, (*itt).first, (*itt).second);

	printf("Old <== New:\n");
	for (Tshiftmap::const_iterator it = new2old.begin(); it != new2old.end(); ++it)
		printf("    [%s, %d] <== %s\n", (*it).second.first, (*it).second.second, (*it).first);
}

void SAtoms::substituteAllLagsAndLeads(FormulaParser& fp, AtomSubstitutions& as)
{
	const char* name;

	int mlead, mlag;
	endovarspan(mlead, mlag);

	// substitute all endo lagged more than 1
	while (NULL != (name = findEndoWithLeadInInterval(mlag, -2)))
		makeAuxVariables(name, -1, -2, mlag, fp, as);
	// substitute all endo leaded more than 1
	while (NULL != (name = findEndoWithLeadInInterval(2, mlead)))
		makeAuxVariables(name, 1, 2, mlead, fp, as);

	exovarspan(mlead, mlag);
	
	// substitute all lagged exo
	while (NULL != (name = findExoWithLeadInInterval(mlag, -1)))
		makeAuxVariables(name, -1, -1, mlag, fp, as);
	// substitute all leaded exo
	while (NULL != (name = findExoWithLeadInInterval(1, mlead)))
		makeAuxVariables(name, 1, 1, mlead, fp, as);

	// notify that substitution have been finished
	as.substitutions_finished(order_type);
}

void SAtoms::substituteAllLagsAndExo1Leads(FormulaParser& fp, AtomSubstitutions& as)
{
	const char* name;

	int mlead, mlag;
	endovarspan(mlead, mlag);

	// substitute all endo lagged more than 1
	while (NULL != (name = findEndoWithLeadInInterval(mlag, -2)))
		makeAuxVariables(name, -1, -2, mlag, fp, as);

	exovarspan(mlead, mlag);
	
	// substitute all lagged exo
	while (NULL != (name = findExoWithLeadInInterval(mlag, -1)))
		makeAuxVariables(name, -1, -1, mlag, fp, as);
	// substitute all leaded exo by 1
	while (NULL != (name = findExoWithLeadInInterval(1,1)))
		makeAuxVariables(name, 1, 1, 1, fp, as);

	// notify that substitution have been finished
	as.substitutions_finished(order_type);
}

const char* SAtoms::findNameWithLeadInInterval(const vector<const char*>& names,
											   int ll1, int ll2) const
{
	for (unsigned int i = 0; i < names.size(); i++) {
		const char* name = names[i];
		DynamicAtoms::Tvarmap::const_iterator it = vars.find(name);
		if (it != vars.end()) {
			const DynamicAtoms::Tlagmap& lmap = (*it).second;
			for (DynamicAtoms::Tlagmap::const_iterator itt = lmap.begin();
				 itt != lmap.end(); ++itt)
				if ((*itt).first >= ll1 && (*itt).first <= ll2)
					return name;
		}
	}

	// nothing found
	return NULL;
}

void SAtoms::attemptAuxName(const char* str, int ll, string& out) const
{
	char c = (ll >= 0)? ((ll == 0)? 'e' : 'p' ) : 'm';
	char absll[100];
	sprintf(absll, "%d", std::abs(ll));
	int iter = 1;
	do {
		out = string(str) + '_';
		for (int i = 0; i < iter; i++)
			out += c;
		if (ll != 0)
			out += absll;
		iter++;
	} while (varnames.query(out.c_str()));
}

void SAtoms::makeAuxVariables(const char* name, int step, int start, int limit_lead,
							  FormulaParser& fp, AtomSubstitutions& as)
{
	if (! (step == 1 || step == -1))
		throw ogu::Exception(__FILE__,__LINE__,
							 "Wrong value of step in SAtoms::makeAuxVariables");
	if (step*start > step*limit_lead)
		throw ogu::Exception(__FILE__,__LINE__,
							 "Wrong value of start in SAtoms::makeAuxVariables");

	// make sure that we do not go further than necessary, this is
	// that the limit lead is not behind maxlead or minlag
	int mlead, mlag;
	varspan(name, mlead, mlag);
	if (step == -1)
		limit_lead = std::max(limit_lead, mlag);
	else
		limit_lead = std::min(limit_lead, mlead);

	// Comment to comments: name="a"; start=-3; step=-1;

	char tmp[500];

	// recover tree index of a previous atom, i.e. set tprev to a tree
	// index of atom "a(-2)"
	int tprev = index(name, start-step);
	if (tprev == -1) {
		sprintf(tmp, "%s(%d)", name, start-step);
		tprev = fp.add_nulary(tmp);
	}

	int ll = start;
	do {
		// either create atom "a_m2(0)" with tree index taux and add
		// equation "a_m2(0)=a(-2)"
		// or 
        // check if "a_m2(0)" has not been already created (with
        // different step), in this case do not add equation "a_m2(0)
        // = a(-2)"
		const char* newname;
		string newname_str;
		int taux;
		if (NULL == (newname=as.get_new4old(name, ll-step))) {
			attemptAuxName(name, ll-step, newname_str);
			newname = newname_str.c_str();
			register_uniq_endo(newname);
			newname = varnames.query(newname);
			sprintf(tmp, "%s(0)", newname);
			taux = fp.add_nulary(tmp);
			// add to substitutions
			as.add_substitution(newname, name, ll-step);

			// add equation "a_m2(0) = a(-2)", this is taux = tprev
			fp.add_formula(fp.add_binary(MINUS, taux, tprev));
		} else {
			// example: exogenous EPS and occurrence at both EPS(-1)
			// EPS(+1)
            // first call makeAuxVariables("EPS",1,1,...) will make endo EPS_p0 = EPS
            // second call makeAuxVariables("EPS",-1,-1,...) will use this EPS_p0
			//             to substitute for EPS(-1)
			taux = index(newname, 0);
			if (taux < 0)
				throw ogu::Exception(__FILE__,__LINE__,
									 "Couldn't find tree index of previously substituted variable");
		}
			
		// create atom "a_m2(-1)" or turn "a(-3)" if any to "a_m2(-1)"; tree index t
		int t = index(name, ll);
		if (t == -1) {
			// no "a(-3)", make t <-> a_m2(-1)
			sprintf(tmp, "%s(%d)", newname, step);
			t = fp.add_nulary(tmp);
		} else {
			// turn a(-3) to a_m2(-1)
			unassign_variable(name, ll, t);
			assign_variable(newname, step, t);
		}

		// next iteration starts with tprev <-> "a_m2(-1)" (this will be made equal to "a_m3(0)")
		tprev = t;
		
		ll += step;
	} while (step*ll <= step*limit_lead);
}
