// Copyright (C) 2006, Ondra Kamenik

// $Id: atom_substitutions.h 42 2007-01-22 21:53:24Z ondra $

#ifndef OGP_ATOM_SUBSTITUTIONS_H
#define OGP_ATOM_SUBSTITUTIONS_H

#include "fine_atoms.h"

#include <string>

namespace ogp {

	using std::string;
	using std::map;
	using std::pair;

	/** This class tracts an information about the performed
	 * substitutions. In fact, there is only one number to keep track
	 * about, this is a number of substitutions. */
	struct SubstInfo {
		int num_substs;
		SubstInfo() : num_substs(0) {}
	};

	/** This class tracks all atom substitutions during the job and
	 * then builds structures when all substitutions are finished. */
	class AtomSubstitutions {
	public:
		typedef pair<const char*, int> Tshiftname;
		typedef map<const char*, Tshiftname, ltstr> Tshiftmap;
		typedef set<Tshiftname> Tshiftnameset;
		typedef map<const char*, Tshiftnameset, ltstr> Toldnamemap;
	protected:
		/** This maps a new name to a shifted old name. This is, one
		 * entry looks as "a_m3 ==> a(-3)", saying that a variable
		 * "a_m3" corresponds to a variable "a" lagged by 3. */
		Tshiftmap new2old;
		/** This is inverse to new2old, which is not unique. For old
		 * name, say "a", it says what new names are derived with what
		 * shifts from the "a". For example, it can map "a" to a two
		 * element set {["a_m3", +3], ["a_p2", -2]}. This says that
		 * leading "a_m3" by 3 one gets old "a" and lagging "a_p2" by
		 * 2 one gets also old "a". */
		Toldnamemap old2new;
		/** This is a reference to old atoms with multiple leads and
		 * lags. They are supposed to be used with parsing finished
		 * being had called, so that the external ordering is
		 * available. */
		const FineAtoms& old_atoms;
		/** This is a reference to new atoms. All name pointers point
		 * to storage of these atoms. */
		FineAtoms& new_atoms;
		/** Substitutions information. */
		SubstInfo info;
	public:
		/** Create the object with reference to the old and new
		 * atoms. In the beginning, old atoms are supposed to be with
		 * parsing_finished() called, and new atoms a simple copy of
		 * old atoms. The new atoms will be an instance of SAtoms. All
		 * substitution job is done by a substitution method of the
		 * new atoms. */
		AtomSubstitutions(const FineAtoms& oa, FineAtoms& na)
			: old_atoms(oa), new_atoms(na) {}
		/** Construct a copy of the object using a different instances
		 * of old atoms and new atoms, which are supposed to be
		 * semantically same as the atoms from as. */
		AtomSubstitutions(const AtomSubstitutions& as, const FineAtoms& oa, FineAtoms& na);
		virtual ~AtomSubstitutions() {}
		/** This is called during the substitution job from the
		 * substitution method of the new atoms. This says that the
		 * new name, say "a_m3" is a substitution of old name "a"
		 * shifted by -3. */
		void add_substitution(const char* newname, const char* oldname, int tshift);
		/** This is called when all substitutions are finished. This
		 * forms the new external ordering of the new atoms and calls
		 * parsing_finished() for the new atoms with the given ordering type. */
		void substitutions_finished(VarOrdering::ord_type ot);
		/** Returns a new name for old name and given tshift. For "a"
		 * and tshift=-3, it returns "a_m3". If there is no such
		 * substitution, it returns NULL. */
		const char* get_new4old(const char* oldname, int tshift) const;
		/** Return new2old. */
		const Tshiftmap& get_new2old() const
			{return new2old;}
		/** Return old2new. */
		const Toldnamemap& get_old2new() const
			{return old2new;}
		/** Return substitution info. */
		const SubstInfo& get_info() const
			{return info;}
		/** Return old atoms. */
		const FineAtoms& get_old_atoms() const
			{return old_atoms;}
		/** Return new atoms. */
		const FineAtoms& get_new_atoms() const
			{return new_atoms;}
		/** Debug print. */
		void print() const;
	};

	class SAtoms : public FineAtoms {
	public:
		SAtoms()
			: FineAtoms() {}
		SAtoms(const SAtoms& sa)
			: FineAtoms(sa) {}
		virtual ~SAtoms() {}
		/** This substitutes all lags and leads for all exogenous and
		 * all lags and leads greater than 1 for all endogenous
		 * variables. This is useful for perfect foresight problems
		 * where we can do that. */
		void substituteAllLagsAndLeads(FormulaParser& fp, AtomSubstitutions& as);
		/** This substitutes all lags of all endo and exo and one step
		 * leads of all exo variables. This is useful for stochastic
		 * models where we cannot solve leads more than 1. */
		void substituteAllLagsAndExo1Leads(FormulaParser& fp, AtomSubstitutions& as);
	protected:
		/** This finds an endogenous variable name which occurs between
		 * ll1 and ll2 included. */
		const char* findEndoWithLeadInInterval(int ll1, int ll2) const
			{return findNameWithLeadInInterval(get_endovars(), ll1, ll2);}
		/** This finds an exogenous variable name which occurs between
		 * ll1 and ll2 included. */
		const char* findExoWithLeadInInterval(int ll1, int ll2) const
			{return findNameWithLeadInInterval(get_exovars(), ll1, ll2);}

		/** This attempts to find a non registered name of the form
		 * <str>_m<abs(ll)> or <str>_p<abs(ll)>. A letter 'p' is
		 * chosen if ll is positive, 'm' if negative. If a name of
		 * such form is already registered, one more character (either
		 * 'p' or 'm') is added and the test is performed again. The
		 * resulting name is returned in a string out. */ 
		void attemptAuxName(const char* str, int ll, string& out) const;

		/** This makes auxiliary variables to eliminate all leads/lags
		 * greater/less than or equal to start up to the limit_lead
		 * for a variable with the given name. If the limit_lead is
		 * greater/less than the maxlead/minlag of the variable, than
		 * maxlead/minlag is used. This process is recorded in
		 * AtomSubstitutions. The new auxiliary variables and their
		 * atoms are created in this object. The auxiliary equations
		 * are created in the given FormulaParser. The value of step
		 * is allowed to be either -1 (lags) or +1 (leads). */
		void makeAuxVariables(const char* name, int step, int start, int limit_lead,
							  FormulaParser& fp, AtomSubstitutions& as);
	private:
		/** This is a worker routine for findEndoWithLeadInInterval
		 * and findExoWithLeadInInterval. */
		const char* findNameWithLeadInInterval(const vector<const char*>& names,
											   int ll1, int ll2) const;

	};

};

#endif

// Local Variables:
// mode:C++
// End:
