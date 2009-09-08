// Copyright (C) 2005, Ondra Kamenik

// $Id: dynamic_atoms.h 2269 2008-11-23 14:33:22Z michel $

#ifndef OGP_DYNAMIC_ATOMS_H
#define OGP_DYNAMIC_ATOMS_H

#include "formula_parser.h"

#include <vector>
#include <map>
#include <set>
#include <string>
#include <cstring>

namespace ogp {
	using std::vector;
	using std::map;
	using std::set;
	using std::string;

	struct ltstr {
		bool operator()(const char* a1, const char* a2) const
			{ return strcmp(a1, a2) < 0; }
	};

	/** Class storing names. We will keep names of variables in
	 * various places, and all these pointers will point to one
	 * storage, which will be responsible for allocation and
	 * deallocation. The main function of the class is to allocate
	 * space for names, and return a pointer of the stored name if
	 * required. */
	class NameStorage {
	protected:
		/** Vector of names allocated, this is the storage. */
		vector<char*> name_store;
		/** Map useful to quickly decide if the name is already
		 * allocated or not. */
		set<const char*, ltstr> name_set;
	public:
		NameStorage() {}
		NameStorage(const NameStorage& stor);
		virtual ~NameStorage();
		/** Query for the name. If the name has been stored, it
		 * returns its address, otherwise 0. */
		const char* query(const char* name) const;
		/** Insert the name if it has not been inserted yet, and
		 * return its new or old allocation. */
		const char* insert(const char* name);
		int num() const
			{return (int)name_store.size();}
		const char* get_name(int i) const
			{return name_store[i];}
		/** Debug print. */
		void print() const;
	};

	class Constants : public AtomValues {
	public:
		/** Type for a map mapping tree indices to double values. */
		typedef map<int,double> Tconstantmap;
		typedef map<int,int> Tintintmap;
	protected:
		/** Map mapping a tree index of a constant to its double value. */ 
		Tconstantmap cmap;
	public:
		Constants() {}
		/** Copy constructor. */
		Constants(const Constants& c)
			: cmap(c.cmap), cinvmap(c.cinvmap) {}
		/** Copy constructor registering the constants in the given
		 * tree. The mapping from old tree indices to new ones is
		 * traced in tmap. */
		Constants(const Constants& c, OperationTree& otree, Tintintmap& tmap)
			{import_constants(c, otree, tmap);}
		/** Import constants registering their tree indices in the
		 * given tree. The mapping form old tree indices to new ones
		 * is traced in tmap. */
		void import_constants(const Constants& c, OperationTree& otree, Tintintmap& tmap);
		/** Implements AtomValues interface. This sets the values to
		 * the evaluation tree EvalTree. */
		void setValues(EvalTree& et) const;
		/** This adds a constant with the given tree index. The
		 * constant must be checked previously and asserted that it
		 * does not exist. */
		void add_constant(int t, double val);
		/** Returns true if the tree index is either an hardwired
		 * constant (initial number OperationTree:num_constants in
		 * OperationTree) or the tree index is a registered constant
		 * by add_constant method. */
		bool is_constant(int t) const;
		double get_constant_value(int t) const;
		/** Return -1 if the given string representation of a constant
		 * is not among the constants (double represenations). If it
		 * is, its tree index is returned. */ 
		int check(const char* str) const;
		/** Debug print. */
		void print() const;
		const Tconstantmap& get_constantmap() const
			{return cmap;}
	private:
		/** Inverse map to Tconstantmap. */
		typedef map<double,int> Tconstantinvmap;
		/** This is an inverse map to cmap. This is only used for fast
		 * queries for the existing double constants in check
		 * method and add_constant. */
		Tconstantinvmap cinvmap;
	};

    /** This class is a parent to Atoms classes which distinguish between
	 * constants (numerical literals), and variables with lags and
	 * leads. This abstraction does not distinguish between a parameter
	 * and a variable without lag or lead. In this sense, everything is a
	 * variable.*/
	class DynamicAtoms : public Atoms, public Constants {
	public:
		/** Definition of a type mapping lags to the indices of the variables. */
		typedef map<int,int> Tlagmap;
	protected:
		/** Definition of a type mapping names of the atoms to Tlagmap. */
		typedef map<const char*, Tlagmap, ltstr> Tvarmap;
		/** Definition of a type mapping indices of variables to the variable names. */
		typedef map<int, const char*> Tindexmap;
		/** This is just a storage for variable names, since all other
		 * instances of a variable name just point to the memory
		 * allocated by this object. */
		NameStorage varnames;
		/** This is the map for variables. Each variable name is
		 * mapped to the Tlagmap, which maps lags/leads to the nulary
		 * term indices in the tree. */
		Tvarmap vars;
		/** This is almost inverse map to the vars. It maps variable
		 * indices to the names. A returned name can be in turn used
		 * as a key in vars. */
		Tindexmap indices;

		/** Number of variables. */
		int nv;
		/** Minimum lag, if there is at least one lag, than this is a negative number. */
		int minlag;
		/** Maximum lead, if there is at least one lead, than this is a positive number. */
		int maxlead;
	public:
		/** Construct empty DynamicAtoms. */
		DynamicAtoms();
		DynamicAtoms(const DynamicAtoms& da);
		virtual ~DynamicAtoms() {}
		/** Check the nulary term identified by its string
		 * representation. The nulary term can be either a constant or
		 * a variable. If constant, -1 is returned so that it could be
		 * assigned regardless if the same constant has already
		 * appeared or not. If variable, then -1 is returned only if
		 * the variable has not been assigned an index, otherwise the
		 * assigned index is returned. */ 
		int check(const char* name) const;
		/** Assign the nulary term identified by its string
		 * representation. This method should be called when check()
		 * returns -1. */
		void assign(const char* name, int t);
		/** Return a number of all variables. */
		int nvar() const
			{ return nv; }
		/** Return the vector of variable indices. */
		vector<int> variables() const;
		/** Return max lead and min lag for a variable given by the
		 * index. If a variable cannot be found, the method retursn
		 * the smallest integer as maxlead and the largest integer as
		 * minlag. */
		void varspan(int t, int& mlead, int& mlag) const;
		/** Return max lead and min lag for a variable given by the
		 * name (without lead, lag). The same is valid if the variable
		 * name cannot be found. */
		void varspan(const char* name, int& mlead, int& mlag) const;
		/** Return max lead and min lag for a vector of variables given by the names. */
		void varspan(const vector<const char*>& names, int& mlead, int& mlag) const;
		/** Return true for all tree indices corresponding to a
		 * variable in the sense of this class. (This is parameters,
		 * exo and endo). Since the semantics of 'variable' will be
		 * changed in subclasses, we use name 'named atom'. These are
		 * all atoms but constants. */
		bool is_named_atom(int t) const;
		/** Return index of the variable described by the variable
		 * name and lag/lead. If it doesn't exist, return -1. */
		int index(const char* name, int ll) const;
		/** Return the lag map for the variable name. */
		const Tlagmap& lagmap(const char* name) const;
		/** Return the variable name for the tree index. It throws an
		 * exception if the tree index t is not a named atom. */
		const char* name(int t) const;
		/** Return the lead/lag for the tree index. It throws an
		 * exception if the tree index t is not a named atom. */
		int lead(int t) const;
		/** Return maximum lead. */
		int get_maxlead() const
			{return maxlead;}
		/** Return minimum lag. */
		int get_minlag() const
			{return minlag;}
		/** Return the name storage to allow querying to other
		 * classes. */
		const NameStorage& get_name_storage() const
			{return varnames;}
		/** Assign the variable with a given lead. The varname must be
		 * from the varnames storage. The method checks if the
		 * variable iwht the given lead/lag is not assigned. If so, an
		 * exception is thrown. */
		void assign_variable(const char* varname, int ll, int t);
		/** Unassign the variable with a given lead and given tree
		 * index. The tree index is only provided as a check. An
		 * exception is thrown if the name, ll, and the tree index t
		 * are not consistent. The method also updates nv, indices,
		 * maxlead and minlag. The varname must be from the varnames
		 * storage. */
		void unassign_variable(const char* varname, int ll, int t);
		/** Debug print. */
		void print() const;
	protected:
		/** Do the check for the variable. A subclass may need to
		 * reimplement this so that it could raise an error if the
		 * variable is not among a given list. */
		virtual int check_variable(const char* name) const;
		/** Assign the constant. */
		void assign_constant(const char* name, int t);
		/** Assign the variable. */
		void assign_variable(const char* name, int t);
		/** The method just updates minlag or/and maxlead. Note that
		 * when assigning variables, the update is done when inserting
		 * to the maps, however, if removing a variable, we need to
		 * call this method. */
		void update_minmaxll();
		/** The method parses the string to recover a variable name
		 * and lag/lead ll. The variable name doesn't contain a lead/lag. */
		virtual void parse_variable(const char* in, string& out, int& ll) const = 0;
	public:
		/** Return true if the str represents a double.*/ 
		static bool is_string_constant(const char* str);
	};


	/** This class is a parent of all orderings of the dynamic atoms
	 * of variables which can appear before t, at t, or after t. It
	 * encapsulates the ordering, and the information about the number
	 * of static (appearing only at time t) predetermined (appearing
	 * before t and possibly at t), both (appearing before t and after
	 * t and possibly at t) and forward looking (appearing after t and
	 * possibly at t).
	 *
	 * The constructor takes a list of variable names. The class also
	 * provides mapping from the ordering of the variables in the list
	 * (outer) to the new ordering (at time t) and back.
	 *
	 * The user of the subclass must call do_ordering() after
	 * initialization.
	 *
	 * The class contains a few preimplemented methods for
	 * ordering. The class is used in this way: Make a subclass, and
	 * implement pure virtual do_ordering() by just plugging a
	 * preimplemented method, or plugging your own implementation. The
	 * method do_ordering() is called by the user after the constructor.
	 */
	class VarOrdering {
	protected:
		/** Number of static variables. */
		int n_stat;
		/** Number of predetermined variables. */
		int n_pred;
		/** Number of both variables. */
		int n_both;
		/** Number of forward looking variables. */
		int n_forw;
		/** This is a set of tree indices corresponding to the
		 * variables at all times as they occur in the formulas. In
		 * fact, since this is used only for derivatives, the ordering
		 * of this vector is only important for ordering of the
		 * derivatives, in other contexts the ordering is not
		 * important, so it is rather a set of indices.*/
		vector<int> der_atoms;
		/** This maps tree index of the variable to the position in
		 * the row of the ordering. One should be careful with making
		 * space in the positions for variables not appearing at time
		 * t. For instance in the pred(t-1), both(t-1), stat(t),
		 * pred(t), both(t), forw(t), both(t+1), forw(t+1) ordering,
		 * the variables x(t-1), y(t-1), x(t+1), z(t-1), z(t), and
		 * z(t+1) having tree indices 6,5,4,3,2,1 will be ordered as
		 * follows: y(t-1), x(t-1), z(t-1), [y(t)], [x(t)], z(t),
		 * x(t+1), where a bracketed expresion means non-existent by
		 * occupying a space. The map thus will look as follows:
		 * {5->0, 6->1, 3->2, 2->5, 3->6}. Note that nothing is mapped
		 * to positions 3 and 4. */ 
		map<int,int> positions;
		/** This maps an ordering of the list of variables in
		 * constructor to the new ordering (at time t). The length is
		 * the number of variables. */
		vector<int> outer2y;
		/** This maps a new ordering to the ordering of the list of
		 * variables in constructor (at time t). The length is the
		 * number of variables. */
		vector<int> y2outer;
		/** This is just a reference for variable names to keep it
		 * from constructor to do_ordering() implementations. */
		const vector<const char*>& varnames;
		/** This is just a reference to atoms to keep it from
		 * constructor to do_ordering() implementations. */
		const DynamicAtoms& atoms;
	public:
		/** This is an enum type for an ordering type implemented by
		 * do_general. */
		enum ord_type {pbspbfbf, bfspbfpb};
		/** Construct the ordering of the variables given by the names
		 * with their dynamic occurrences defined by the atoms. It
		 * calls the virtual method do_ordering which can be
		 * reimplemented. */
		VarOrdering(const vector<const char*>& vnames, const DynamicAtoms& a)
			: n_stat(0), n_pred(0), n_both(0), n_forw(0), varnames(vnames), atoms(a)
			{}
		VarOrdering(const VarOrdering& vo, const vector<const char*>& vnames,
					const DynamicAtoms& a);
		virtual VarOrdering* clone(const vector<const char*>& vnames,
								   const DynamicAtoms& a) const = 0;
		/** Destructor does nothing here. */
		virtual ~VarOrdering() {}
		/** This is the method setting the ordering and the map. A
		 * subclass must reimplement it, possibly using a
		 * preimplemented ordering. This method must be called by the
		 * user after the class has been created. */
		virtual void do_ordering() = 0;
		/** Return number of static. */
		int nstat() const
			{return n_stat;}
		/** Return number of predetermined. */
		int npred() const
			{return n_pred;}
		/** Return number of both. */
		int nboth() const
			{return n_both;}
		/** Return number of forward looking. */
		int nforw() const
			{return n_forw;}
		/** Return the set of tree indices for derivatives. */
		const vector<int>& get_der_atoms() const
			{return der_atoms;}
		/** Return the y2outer. */
		const vector<int>& get_y2outer() const
			{return y2outer;}
		/** Return the outer2y. */
		const vector<int>& get_outer2y() const
			{return outer2y;}
		/** Query the atom given by the tree index. True is returned
		 * if the atom is one of the variables in the object. */
		bool check(int t) const;
		/** Return the position of the atom (nulary term) given by a
		 * tree index. It is a lookup to the map. If the atom cannot
		 * be found, the exception is raised. */
		int get_pos_of(int t) const;
		/** This returns a length of ordered row of atoms. In all
		 * cases so far, it does not depend on the ordering and it is
		 * as follows. */
		int length() const
			{return n_stat+2*n_pred+3*n_both+2*n_forw;}
		/** Debug print. */
		void print() const;
	protected:
		/** This is a general ordering method which orders the
		 * variables by the given ordering ord_type. See documentation
		 * for respective do_ methods. */
		void do_general(ord_type ordering);
		/** This is a preimplemented ordering for do_ordering()
		 * method. It assumes that the variables appear only at time
		 * t-1, t, t+1. It orders the atoms as pred(t-1), both(t-1),
		 * stat(t), pred(t), both(t), forw(t), both(t+1),
		 * forw(t+1). It builds the der_atoms, the map of positions,
		 * as well as y2outer and outer2y. */
		void do_pbspbfbf()
			{do_general(pbspbfbf);}
		/** This is a preimplemented ordering for do_ordering()
		 * method. It assumes that the variables appear only at time
		 * t-1, t, t+1. It orders the atoms as both(t+1), forw(t+1),
		 * stat(t), pred(t), both(t), forw(t), pred(t-1),
		 * both(t-1). It builds the der_atoms, the map of positions,
		 * as well as y2outer and outer2y. */
		void do_bfspbfpb()
			{do_general(bfspbfpb);}
		/** This is a preimplemented ordering for do_ordering()
		 * method. It makes no assumptions about occurences of
		 * variables at different times. It orders the atoms with
		 * increasing time keeping the given ordering within one
		 * time. This implies that y2outer and outer2y will be
		 * identities. The der_atoms will be just a sequence of atoms
		 * from the least to the most time preserving the order of atoms
		 * within one time. */
		void do_increasing_time();
	private:
		/** Declare this copy constructor as private to hide it. */
		VarOrdering(const VarOrdering& vo);
	};

};

#endif

// Local Variables:
// mode:C++
// End:
