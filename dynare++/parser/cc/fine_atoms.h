// Copyright (C) 2005, Ondra Kamenik

// $Id: fine_atoms.h 1759 2008-03-31 14:25:20Z kamenik $

#ifndef OGP_FINE_ATOMS_H
#define OGP_FINE_ATOMS_H

#include "dynamic_atoms.h"

#include <vector>
#include <string>

namespace ogp {

	using std::vector;
	using std::string;

	/** This is just ordering used for endogenous variables. It
	 * assumes that we have only time t-1, t, and t+1, orders them as
	 * pred(t-1), both(t-1), stat(t), pred(t), both(t), forw(t),
	 * both(t+1), forw(t+1). */
	class EndoVarOrdering1 : public VarOrdering {
	public:
		EndoVarOrdering1(const vector<const char*>& vnames, const DynamicAtoms& a)
			: VarOrdering(vnames, a) {}
		EndoVarOrdering1(const EndoVarOrdering1& vo, const vector<const char*>& vnames,
						 const DynamicAtoms& a)
			: VarOrdering(vo, vnames, a) {}
		VarOrdering* clone(const vector<const char*>& vnames, const DynamicAtoms& a) const
			{return new EndoVarOrdering1(*this, vnames, a);}
		void do_ordering()
			{do_pbspbfbf();}
	};

	/** This is just another ordering used for endogenous
	 * variables. It assumes that we have only time t-1, t, and t+1,
	 * orders them as both(t+1), forw(t+1), pred(t-1), both(t-1),
	 * stat(t), pred(t), both(t), forw(t). */
	class EndoVarOrdering2 : public VarOrdering {
	public:
		EndoVarOrdering2(const vector<const char*>& vnames, const DynamicAtoms& a)
			: VarOrdering(vnames, a) {}
		EndoVarOrdering2(const EndoVarOrdering2& vo, const vector<const char*>& vnames,
						 const DynamicAtoms& a)
			: VarOrdering(vo, vnames, a) {}
		VarOrdering* clone(const vector<const char*>& vnames, const DynamicAtoms& a) const
			{return new EndoVarOrdering2(*this, vnames, a);}
		void do_ordering()
			{do_bfspbfpb();}
	};

	/** This is just ordering used for exogenous variables. It makes
	 * no assumptions about their timing. It orders them from the
	 * least time to the latest time. */
	class ExoVarOrdering : public VarOrdering {
	public:
		ExoVarOrdering(const vector<const char*>& vnames, const DynamicAtoms& a)
			: VarOrdering(vnames, a) {}
		ExoVarOrdering(const ExoVarOrdering& vo, const vector<const char*>& vnames,
					   const DynamicAtoms& a)
			: VarOrdering(vo, vnames, a) {}
		VarOrdering* clone(const vector<const char*>& vnames, const DynamicAtoms& a) const
			{return new ExoVarOrdering(*this, vnames, a);}
		void do_ordering()
			{do_increasing_time();}
	};

	class FineAtoms;

	/** This class provides an outer ordering of all variables (endo
	 * and exo). It maps the ordering to the particular outer
	 * orderings of endo and exo. It works tightly with the FineAtoms
	 * class. */
	class AllvarOuterOrdering {
	protected:
		/** Type for a map mapping a variable name to an integer. */
		typedef map<const char*, int, ltstr> Tvarintmap;
		/** Reference to atoms. */
		const FineAtoms& atoms;
		/** The vector of all endo and exo variables in outer
		 * ordering. The pointers point to storage in atoms. */
		vector<const char*> allvar;
		/** The mapping from outer endogenous to outer all. For
		 * example endo2all[0] is the order of the first outer
		 * endogenous variable in the allvar ordering. */
		vector<int> endo2all;
		/** The mapping from outer exogenous to outer all. For example
		 * exo2all[0] is the order of the first outer exogenous
		 * variables in the allvar ordering. */
		vector<int> exo2all;
	public:
		/** Construct the allvar outer ordering from the provided
		 * sequence of endo and exo names. The names can have an
		 * arbitrary storage, the storage is transformed to the atoms
		 * storage. An exception is thrown if either the list is not
		 * exhaustive, or some string is not a variable. */
		AllvarOuterOrdering(const vector<const char*>& allvar_outer, const FineAtoms& a);
		/** Copy constructor using the storage of provided atoms. */
		AllvarOuterOrdering(const AllvarOuterOrdering& allvar_outer, const FineAtoms& a);
		/** Return endo2all mapping. */
		const vector<int>& get_endo2all() const
			{return endo2all;}
		/** Return exo2all mapping. */
		const vector<int>& get_exo2all() const
			{return exo2all;}
		/** Return the allvar ordering. */
		const vector<const char*>& get_allvar() const
			{return allvar;}
	};

	/** This class refines the DynamicAtoms by distinguishing among
	 * parameters (no lag and leads) and endogenous and exogenous
	 * variables (with lags and leads). For parameters, endogenous and
	 * exogenous, it defines outer orderings and internal
	 * orderings. The internal orderings are created by
	 * parsing_finished() method when it is sure that no new variables
	 * would be registered. The outer orderings are given by the order
	 * of calls of registering methods.
     * 
     * In addition, the class also defines outer ordering of
     * endogenous and exogenous variables. This is input as a
     * parameter to parsing_finished(). By default, this whole outer
     * ordering is just a concatenation of outer ordering of
     * endogenous and exogenous variables.
	 *
	 * The internal ordering of all endo and exo variables is just a
	 * concatenation of endo and exo variables in their internal
	 * orderings. This is the ordering with respect to which all
	 * derivatives are taken. */
	class FineAtoms : public DynamicAtoms {
		friend class AllvarOuterOrdering;
	protected:
		typedef map<const char*, int, ltstr> Tvarintmap;
	private:
		/** The vector of parameters names. The order gives the order
		 * the data is communicated with outside world. */
		vector<const char*> params;
		/** A map mapping a name of a parameter to an index in the outer
		 * ordering. */
		Tvarintmap param_outer_map;
		/** The vector of endogenous variables. This defines the order
		 * like parameters. */
		vector<const char*> endovars;
		/** A map mapping a name of an endogenous variable to an index
		 * in the outer ordering. */
		Tvarintmap endo_outer_map;
		/** The vector of exogenous variables. Also defines the order
		 * like parameters and endovars. */
		vector<const char*> exovars;
		/** A map mapping a name of an exogenous variable to an index
		 * in the outer ordering. */
		Tvarintmap exo_outer_map;

	protected:
		/** This is the internal ordering of all atoms corresponding
		 * to endogenous variables. It is constructed by
		 * parsing_finished() method, which should be called after all
		 * parsing jobs have been finished. */ 
		VarOrdering* endo_order;
		/** This is the internal ordering of all atoms corresponding
		 * to exogenous variables. It has the same handling as
		 * endo_order. */
		VarOrdering* exo_order;
		/** This is the all variables outer ordering. It is
		 * constructed by parsing finished. */
		AllvarOuterOrdering* allvar_order;
		/** This vector defines a set of atoms as tree indices used
		 * for differentiation. The order of the atoms in this vector
		 * defines ordering of the derivative tensors. The ordering is
		 * a concatenation of atoms from endo_order and then
		 * exo_order. This vector is setup by parsing_finished() and
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
		FineAtoms()
			: endo_order(NULL), exo_order(NULL), allvar_order(NULL) {}
		FineAtoms(const FineAtoms& fa);
		/** Deletes endo_order and exo_order. */
		virtual ~FineAtoms()
			{
				if (endo_order) delete endo_order;
				if (exo_order) delete exo_order;
				if (allvar_order) delete allvar_order;
			}
		/** Overrides DynamicAtoms::check_variable so that the error
		 * would be raised if the variable name is not declared. A
		 * variable is declared by inserting it to
		 * DynamicAtoms::varnames. This is a responsibility of a
		 * subclass. */
		int check_variable(const char* name) const;
		/** This calculates min lag and max lead of endogenous variables. */
		void endovarspan(int& mlead, int& mlag) const
			{varspan(endovars, mlead, mlag);}
		/** This calculates mim lag and max lead of exogenous variables. */
		void exovarspan(int& mlead, int& mlag) const
			{varspan(exovars, mlead, mlag);}
		/** This calculates the number of periods in which at least
		 * one exogenous variable occurs. */
		int num_exo_periods() const;
		/** Return an (external) ordering of parameters. */
		const vector<const char*>& get_params() const
			{return params;}
		/** Return an external ordering of endogenous variables. */
		const vector<const char*>& get_endovars() const
			{return endovars;}
		/** Return an external ordering of exogenous variables. */
		const vector<const char*>& get_exovars() const
			{return exovars;}
		/** This constructs internal orderings and makes the indices
		 * returned by variables method available. Further it
		 * constructs outer ordering of all variables by a simple
		 * concatenation of outer endogenous and outer exogenous. In
		 * addition, it makes nstat, npred, nboth, nforw available. */
		void parsing_finished(VarOrdering::ord_type ot);
		/** This does the same thing as
		 * parsing_finished(VarOrdering::ord_type) plus it allows for
		 * inputing a different outer ordering of all variables. The
		 * ordering is input as a list of strings, their storage can
		 * be arbitrary. */
		void parsing_finished(VarOrdering::ord_type ot, const vector<const char*> avo);
		/** Return the external ordering of all variables (endo and
		 * exo). This is either the second argument to
		 * parsing_finished or the default external ordering. This
		 * must be called only after parsing_finished. */
		const vector<const char*>& get_allvar() const;
		/** Return the map from outer ordering of endo variables to
		 * the allvar ordering. This must be called only after
		 * parsing_finished. */
		const vector<int>& outer_endo2all() const;
		/** Return the map from outer ordering of exo variables to
		 * the allvar ordering. This must be called only after
		 * parsing_finished. */
		const vector<int>& outer_exo2all() const;
		/** Return the atoms with respect to which we are going to
		 * differentiate. This must be called after
		 * parsing_finished. */
		vector<int> variables() const;
		/** Return the number of static. */
		int nstat() const;
		/** Return the number of predetermined. */
		int npred() const;
		/** Return the number of both. */
		int nboth() const;
		/** Return the number of forward looking. */
		int nforw() const;
		/** Return the index of an endogenous atom given by tree index in
		 * the endo ordering. This must be also called only after
		 * parsing_finished(). */
		int get_pos_of_endo(int t) const;
		/** Return the index of an exogenous atom given by tree index in
		 * the exo ordering. This must be also called only after
		 * parsing_finished(). */
		int get_pos_of_exo(int t) const;
		/** Return the index of either endogenous or exogenous atom
		 * given by tree index in the concatenated ordering of
		 * endogenous and exogenous atoms. This must be also called
		 * only after parsing_finished(). */
		int get_pos_of_all(int t) const;
		/** Return the mapping from endogenous at time t to outer
		 * ordering of endogenous. */
		const vector<int>& y2outer_endo() const;
		/** Return the mapping from the outer ordering of endogenous to endogenous
		 * at time t. */
		const vector<int>& outer2y_endo() const;
		/** Return the mapping from exogenous at time t to outer
		 * ordering of exogenous. */
		const vector<int>& y2outer_exo() const;
		/** Return the mapping from the outer ordering of exogenous to exogenous
		 * at time t. */
		const vector<int>& outer2y_exo() const;
		/** Return the endo_atoms_map. */
		const vector<int>& get_endo_atoms_map() const;
		/** Return the exo_atoms_map. */
		const vector<int>& get_exo_atoms_map() const;
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
		/** Return an index in the outer ordering of all variables
		 * (endo and exo) for a given name. An exception is thrown if
		 * the name is not a variable. This must be called only after
		 * parsing_finished(). */
		int name2outer_allvar(const char* name) const;
		/** Return the number of endogenous variables at time t-1, these are state
		 * variables. */
		int nys() const
			{return npred()+nboth();}
		/** Return the number of endogenous variables at time t+1. */
		int nyss() const
			{return nboth()+nforw();}
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
		/** This performs the common part of parsing_finished(), which
		 * is a construction of internal orderings. */
		void make_internal_orderings(VarOrdering::ord_type ot);
	protected:
		/** This remembers the ordering type of the last call make_internal_ordering. */
		VarOrdering::ord_type order_type;
	};
};

#endif

// Local Variables:
// mode:C++
// End:
