// Copyright (C) 2006, Ondra Kamenik

// $Id$

#ifndef FORW_SUBST_BUILDER_H
#define FORW_SUBST_BUILDER_H


#include "dynare_model.h"

namespace ogdyn {

	/** This struct encapsulates information about the process of
	 * forward substitutions. */
	struct ForwSubstInfo {
		int num_affected_equations;
		int num_subst_terms;
		int num_aux_variables;
		int num_new_terms;
		ForwSubstInfo()
			: num_affected_equations(0),
			  num_subst_terms(0),
			  num_aux_variables(0),
			  num_new_terms(0) {}
	};

	class ForwSubstBuilder {
		typedef map<int, const char*> Ttermauxmap;
	protected:
		/** Reference to the model, to which we will add equations and
		 * change some equations. */
		DynareModel& model;
		/** A map mapping new auxiliary variables to the terms in the
		 * tree in the DynareModel. */
		Tsubstmap aux_map;
		/** Information about the substitutions. */
		ForwSubstInfo info;
	public:
		/** Do all the jobs needed. This scans all equations in the
		 * model, and for equations containing forward looking
		 * variables greater than 1 lead, it makes corresponding
		 * substitutions. Basically, it breaks each equation to its
		 * non-linear components and creates substitutions for these
		 * components, not for whole equation. This is because the
		 * expectation operator can go through the linear part of the
		 * function. This will save us many occurrences of other
		 * variables involved in the equation. */
		ForwSubstBuilder(DynareModel& m);
		/** Copy constructor with a new instance of the model. */
		ForwSubstBuilder(const ForwSubstBuilder& b, DynareModel& m);
		/** Return the auxiliary variable mapping. */
		const Tsubstmap& get_aux_map() const
			{return aux_map;}
		/** Return the information. */
		const ForwSubstInfo& get_info() const
			{return info;}
	private:
		ForwSubstBuilder(const ForwSubstBuilder& b);
		/** This method takes a nonlinear term t, and if it has leads
		 * of greater than 1, then it substitutes the term for the new
		 * variable (or string of variables). Note that the
		 * substitution is done by DynamicAtoms::assign_variable. This
		 * means that the substitution is made for all other
		 * ocurrences of t in the model. So there is no need of
		 * tracking already substituted terms. The other two
		 * parameters are just for identification of the new auxiliary
		 * variables. When called from the constructor, i is an
		 * equation number, j is an order of the non-linear term in
		 * the equation. */
		void substitute_for_term(int t, int i, int j);
		/** This is called just at the end of the job. It unassigns
		 * all nulary terms with a lead greater than 1. */
		void unassign_gt_1_leads();
		/** This unassigns all leads greater than 1 of the given name. */
		void unassign_gt_1_leads(const char* name);
	};
};

#endif

// Local Variables:
// mode:C++
// End:
