// Copyright (C) 2006, Ondra Kamenik

// $Id: static_fine_atoms.h 42 2007-01-22 21:53:24Z ondra $

#ifndef OGP_STATIC_FINE_ATOMS_H
#define OGP_STATIC_FINE_ATOMS_H

#include "static_atoms.h"
#include "fine_atoms.h"

namespace ogp {

	/** This class represents static atoms distinguishing between
	 * parameters, endogenous and exogenous variables. The class
	 * maintains also ordering of all three categories (referenced as
	 * outer or inner, since there is only one ordering). It can be
	 * constructed either from scratch, or from fine dynamic atoms. In
	 * the latter case, one can decide if the ordering of this static
	 * atoms should be internal or external ordering of the original
	 * dynamic fine atoms. */
	class StaticFineAtoms : public StaticAtoms {
	public:
		typedef map<int,int> Tintintmap;
	protected:
		typedef map<const char*, int, ltstr> Tvarintmap;
	private:
		/** The vector of parameter names, gives the parameter
		 * ordering. */
		vector<const char*> params;
		/** A map mappping a parameter name to an index in the ordering. */
		Tvarintmap param_outer_map;
		/** The vector of endogenous variables. This defines the order
		 * like parameters. */
		vector<const char*> endovars;
		/** A map mapping a name of an endogenous variable to an index
		 * in the ordering. */
		Tvarintmap endo_outer_map;
		/** The vector of exogenous variables. Also defines the order
		 * like parameters and endovars. */
		vector<const char*> exovars;
		/** A map mapping a name of an exogenous variable to an index
		 * in the outer ordering. */
		Tvarintmap exo_outer_map;
		/** This vector defines a set of atoms as tree indices used
		 * for differentiation. The order of the atoms in is the
		 * concatenation of the outer ordering of endogenous and
		 * exogenous. This vector is setup by parsing_finished() and
		 * is returned by variables(). */
		vector<int> der_atoms;
		/** This is a mapping from endogenous atoms to all atoms in
		 * der_atoms member. The mapping maps index in endogenous atom
		 * ordering to index (not value) in der_atoms. It is useful if
		 * one wants to evaluate derivatives wrt only endogenous
		 * variables. It is set by parsing_finished(). By definition,
		 * it is monotone. */
		vector<int> endo_atoms_map;
		/** This is a mapping from exogenous atoms to all atoms in
		 * der_atoms member. It is the same as endo_atoms_map for
		 * atoms of exogenous variables. */
		vector<int> exo_atoms_map;		
	public:
		StaticFineAtoms() {}
		/** Copy constructor making a new storage for atom names. */
		StaticFineAtoms(const StaticFineAtoms& sfa);
		/** Conversion from dynamic FineAtoms taking its outer
		 * ordering as ordering of parameters, endogenous and
		 * exogenous. A biproduct is an integer to integer map mapping
		 * tree indices of the dynamic atoms to tree indices of the
		 * static atoms. */
		StaticFineAtoms(const FineAtoms& fa, OperationTree& otree, Tintintmap& tmap)
			{StaticFineAtoms::import_atoms(fa, otree, tmap);}
		/** Conversion from dynamic FineAtoms taking its internal
		 * ordering as ordering of parameters, endogenous and
		 * exogenous. A biproduct is an integer to integer map mapping
		 * tree indices of the dynamic atoms to tree indices of the
		 * static atoms. */
		StaticFineAtoms(const FineAtoms& fa, OperationTree& otree, Tintintmap& tmap,
						const char* dummy)
			{StaticFineAtoms::import_atoms(fa, otree, tmap, dummy);}
		virtual ~StaticFineAtoms() {}
		/** This adds atoms from dynamic atoms inserting new tree
		 * indices to the given tree and tracing the mapping from old
		 * atoms to new atoms in tmap. The ordering of the static
		 * atoms is the same as outer ordering of dynamic atoms. */
		void import_atoms(const FineAtoms& fa, OperationTree& otree, Tintintmap& tmap);
		/** This adds atoms from dynamic atoms inserting new tree
		 * indices to the given tree and tracing the mapping from old
		 * atoms to new atoms in tmap. The ordering of the static
		 * atoms is the same as internal ordering of dynamic atoms. */
		void import_atoms(const FineAtoms& fa, OperationTree& otree, Tintintmap& tmap,
						  const char* dummy);
		/** Overrides StaticAtoms::check_variable so that the error
		 * would be raised if the variable name is not declared. A
		 * variable is declared by inserting it to
		 * StaticAtoms::varnames, which is done with registering
		 * methods. This a responsibility of a subclass. */
		int check_variable(const char* name) const;
		/** Return an (external) ordering of parameters. */
		const vector<const char*>& get_params() const
			{return params;}
		/** Return an external ordering of endogenous variables. */
		const vector<const char*>& get_endovars() const
			{return endovars;}
		/** Return an external ordering of exogenous variables. */
		const vector<const char*>& get_exovars() const
			{return exovars;}
		/** This constructs der_atoms, and the endo_endoms_map and
		 * exo_atoms_map, which can be created only after the parsing
		 * is finished. */
		void parsing_finished();
		/** Return the atoms with respect to which we are going to
		 * differentiate. */
		vector<int> variables() const
			{return der_atoms;}
		/** Return the endo_atoms_map. */
		const vector<int>& get_endo_atoms_map() const
			{return endo_atoms_map;}
		/** Return the exo_atoms_map. */
		const vector<int>& get_exo_atoms_map() const
			{return endo_atoms_map;}
		/** Return an index in the outer ordering of a given
		 * parameter. An exception is thrown if the name is not a
		 * parameter. */
		int name2outer_param(const char* name) const;
		/** Return an index in the outer ordering of a given
		 * endogenous variable. An exception is thrown if the name is not a
		 * and endogenous variable. */
		int name2outer_endo(const char* name) const;
		/** Return an index in the outer ordering of a given
		 * exogenous variable. An exception is thrown if the name is not a
		 * and exogenous variable. */
		int name2outer_exo(const char* name) const;
		/** Return the number of endogenous variables. */
		int ny() const
			{return endovars.size();}
		/** Return the number of exogenous variables. */
		int nexo() const
			{return (int)exovars.size();}
		/** Return the number of parameters. */
		int np() const
			{return (int)(params.size());}
		/** Register unique endogenous variable name. The order of
		 * calls defines the endo outer ordering. The method is
		 * virtual, since a superclass may want to do some additional
		 * action. */
		virtual void register_uniq_endo(const char* name);
		/** Register unique exogenous variable name. The order of
		 * calls defines the exo outer ordering. The method is
		 * virtual, since a superclass may want to do somem additional
		 * action. */
		virtual void register_uniq_exo(const char* name);
		/** Register unique parameter name. The order of calls defines
		 * the param outer ordering. The method is
		 * virtual, since a superclass may want to do somem additional
		 * action. */
		virtual void register_uniq_param(const char* name);
		/** Debug print. */
		void print() const;
	private:
		/** Add endogenous variable name, which is already in the name
		 * storage. */
		void register_endo(const char* name);
		/** Add exogenous variable name, which is already in the name
		 * storage. */
		void register_exo(const char* name);
		/** Add parameter name, which is already in the name
		 * storage. */
		void register_param(const char* name);
	};

};

#endif

// Local Variables:
// mode:C++
// End:
