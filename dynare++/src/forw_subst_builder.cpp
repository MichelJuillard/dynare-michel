// Copyright (C) 2006-2011, Ondra Kamenik

#include "forw_subst_builder.h"

using namespace ogdyn;

ForwSubstBuilder::ForwSubstBuilder(DynareModel& m)
	: model(m)
{
	info.num_new_terms -= model.getParser().getTree().get_num_op();

	// go through all equations
	int neq = model.eqs.nformulas();
	for (int i = 0; i < neq; i++) {
		int ft = model.eqs.formula(i);
		int mlead, mlag;
		model.termspan(ft, mlead, mlag);
		// if equation is too forward looking
		if (mlead > 1) {
			info.num_affected_equations++;
			// break it to non-linear terms
			unordered_set<int> nlt = model.get_nonlinear_subterms(ft);
			int j = 0; // indexes subterms
			// and make substitutions for all these non-linear subterms
			for (unordered_set<int>::const_iterator it = nlt.begin();
				 it != nlt.end(); ++it, ++j)
				substitute_for_term(*it, i, j);
		}
	}
	// unassign all variables with lead greater than 1
	unassign_gt_1_leads();

	// forget the derivatives in the tree because some variables could
	// have been unassigned
	model.eqs.getTree().forget_derivative_maps();

	info.num_new_terms += model.getParser().getTree().get_num_op();
}

void ForwSubstBuilder::substitute_for_term(int t, int i, int j)
{
	int mlead, mlag;
	model.termspan(t, mlead, mlag);
	if (mlead > 1) {
		info.num_subst_terms++;
		// Example for comments: let t = f(x(+4))
		// first make lagsubst be substitution setting f(x(+4)) to f(x(+1))
		// this is lag = -3 (1-mlead)
		map<int,int> lagsubst;
		model.variable_shift_map(model.eqs.nulary_of_term(t), 1-mlead, lagsubst);
		int lagt = model.eqs.add_substitution(t, lagsubst);
		// now maxlead of lagt is +1
		// add AUXLD_*_*_1 = f(x(+1)) to the model
		char name[100];
		sprintf(name, "AUXLD_%d_%d_%d", i, j, 1);
		model.atoms.register_uniq_endo(name);
		info.num_aux_variables++;
		const char* ss = model.atoms.get_name_storage().query(name);
		int auxt = model.eqs.add_nulary(name);
		model.eqs.add_formula(model.eqs.add_binary(ogp::MINUS, auxt, lagt));
		aux_map.insert(Tsubstmap::value_type(ss, lagt));
		// now add variables and equations
		// AUXLD_*_*_2 = AUXLD_*_*_1(+1) through
		// AUXLD_*_*_{mlead-1} = AUXLD_*_*_{mlead-2}(+1)
		for (int ll = 1; ll <= mlead-2; ll++) {
			// create AUXLD_*_*_{ll}(+1)
			sprintf(name, "AUXLD_%d_%d_%d(+1)", i, j, ll);
			int lastauxt_lead = model.eqs.add_nulary(name);
			// create AUXLD_*_*{ll+1}
			sprintf(name, "AUXLD_%d_%d_%d", i, j, ll+1);
			model.atoms.register_uniq_endo(name);
			info.num_aux_variables++;
			ss = model.atoms.get_name_storage().query(name);
			auxt = model.eqs.add_nulary(name);
			// add AUXLD_*_*_{ll+1} = AUXLD_*_*_{ll}(+1)
			model.eqs.add_formula(model.eqs.add_binary(ogp::MINUS, auxt, lastauxt_lead));
			// add substitution to the map; todo: this
			// works well because in the context where
			// aux_map is used the timing doesn't matter,
			// however, it is misleading, needs to be
			// changed
			aux_map.insert(Tsubstmap::value_type(ss, lagt));
		}

		// now we have to substitute AUXLEAD_*_*{mlead-1}(+1) for t
		model.substitute_atom_for_term(ss, +1, t);
	}
}

void ForwSubstBuilder::unassign_gt_1_leads(const char* name)
{
	const char* ss = model.atoms.get_name_storage().query(name);
	int mlead, mlag;
	model.atoms.varspan(name, mlead, mlag);
	for (int ll = 2; ll <= mlead; ll++) {
		int t = model.atoms.index(ss, ll);
		if (t != -1)
			model.atoms.unassign_variable(ss, ll, t);
	}
}

void ForwSubstBuilder::unassign_gt_1_leads()
{
	const vector<const char*>& endovars = model.atoms.get_endovars();
	for (unsigned int i = 0; i < endovars.size(); i++)
		unassign_gt_1_leads(endovars[i]);
	const vector<const char*>& exovars = model.atoms.get_exovars();
	for (unsigned int i = 0; i < exovars.size(); i++)
		unassign_gt_1_leads(exovars[i]);	
}

ForwSubstBuilder::ForwSubstBuilder(const ForwSubstBuilder& b, DynareModel& m)
	: model(m)
{
	for (Tsubstmap::const_iterator it = b.aux_map.begin();
		 it != b.aux_map.end(); ++it) {
		const char* ss = m.atoms.get_name_storage().query((*it).first);
		aux_map.insert(Tsubstmap::value_type(ss, (*it).second));
	}
}
