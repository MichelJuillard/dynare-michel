// Copyright (C) 2005, Ondra Kamenik

// $Id: formula_parser.h 1760 2008-03-31 14:26:35Z kamenik $

#ifndef OGP_FORMULA_PARSER_H
#define OGP_FORMULA_PARSER_H

#include "tree.h"

namespace ogp {
	using std::vector;

	/** Pure virtual class defining a minimal interface for
	 * representation of nulary terms within FormulaParser. */
	class Atoms {
	public:
		Atoms() {}
		virtual ~Atoms() {}
		/** This returns previously assigned internal index to the
		 * given atom, or returns -1 if the atom has not been assigned
		 * yet. The method can raise an exception, if the Atoms
		 * implementation is strict and the name is not among
		 * prescribed possible values. */
		virtual int check(const char* name) const = 0;
		/** This method assigns an internal index to the nulary term
		 * described by the name. The internal index is allocated by
		 * OperationTree class. */
		virtual void assign(const char* name, int t) = 0;
		/** Returns a number of variables which will be used for
		 * differentiations. */
		virtual int nvar() const = 0;
		/** Returns a vector of variable's internal indices which will
		 * be used for differentiations. */
		virtual vector<int> variables() const = 0;
		/** Debug print. */
		virtual void print() const = 0;
	};

	/** Pure virtual class defining interface for all classes able to
	 * set nulary terms to evaluation tree EvalTree. The
	 * implementations of this class will have to be connected with
	 * Atoms to have knowledge about the atoms and their indices in
	 * the tree, and will call EvalTree::set_nulary. */
	class AtomValues {
	public:
		virtual ~AtomValues() {}
		virtual void setValues(EvalTree& et) const = 0;
	};

	class FormulaDerEvaluator;
	class FoldMultiIndex;
	/** For ordering FoldMultiIndex in the std::map. */
	struct ltfmi {
		bool operator()(const FoldMultiIndex& i1, const FoldMultiIndex& i2) const;
	};

	/** This class stores derivatives (tree indices) of one formula
	 * for all orders upto a given one. It stores the derivatives as a
	 * sequence (vector) of these tree indices and sequence of the
	 * multidimensional indices of variables wrt which the derivatives
	 * were taken. In order to speed up querying for a derivative
	 * given the variables, we have a map mapping the multidimensional
	 * index to the order of the derivative in the sequence.
	 * 
	 * The only reason we do not have only this map is that the
	 * iterators of the map do not survive the insertions to the map,
	 * and implementation of the constructor has to be very difficult.
	 */
	class FormulaDerivatives {
		friend class FormulaDerEvaluator;
	protected:
		/** Vector of derivatives. This is a list of derivatives (tree
		 * indices), the ordering is given by the algorithm used to
		 * create it. Currently, it starts with zero-th derivative,
		 * the formula itself and carries with first order, second,
		 * etc. */
		vector<int> tder;
		/** Vector of multiindices corresponding to the vector of
		 * derivatives. */
		vector<FoldMultiIndex> indices;
		/** For retrieving derivatives via a multiindex, we have a map
		 * mapping a multiindex to a derivative in the tder
		 * ordering. This means that indices[ind2der[index]] == index. */
		typedef map<FoldMultiIndex, int, ltfmi> Tfmiintmap;
		Tfmiintmap ind2der;
		/** The number of variables. */
		int nvar;
		/** The maximum order of derivatives. */
		int order;
	public:
		/** The constructor allocates and fills the sequence of the
		 * indices of derivatives for a formula.
		 * @param otree the OperationTree for which all work is done
		 * and to which the derivatives are added.
		 * @param vars the vector of nulary terms in the tree; the
		 * derivatives are taken with respect to these variables in
		 * the ordering given by the vector.
		 * @param f the index of the formula being differentiated. The
		 * zero derivative is set to f.
		 * @param max_order the maximum order of differentiation.
		 */ 
		FormulaDerivatives(OperationTree& otree, const vector<int>& vars, int f, int max_order);
		/** Copy constructor. */
		FormulaDerivatives(const FormulaDerivatives& fd);
		virtual ~FormulaDerivatives(){}
		/** Random access to the derivatives via multiindex. */
		int derivative(const FoldMultiIndex& mi) const;
		/** Return the order. */
		int get_order() const
			{return order;}
		/** Debug print. */
		void print(const OperationTree& otree) const;
	};

	class FormulaEvaluator;

	/** This class is able to parse a number of formulas and
	 * differentiate them. The life cycle of the object is as follows:
	 * After it is created, a few calls to parse will add formulas
	 * (zero derivatives) to the object. Then a method differentiate()
	 * can be called and a vector of pointers to derivatives for each
	 * formula is created. After this, no one should call other
	 * parse() or differentiate(). A const reference of the object can
	 * be used in constructors of FormulaEvaluator and
	 * FormulaDerEvaluator in order to evaluate formulas (zero
	 * derivatives) and higher derivatives resp. */
	class FormulaParser {
		friend class FormulaCustomEvaluator;
		friend class FormulaDerEvaluator;
	protected:
		/** The OperationTree of all formulas, including derivatives. */
		OperationTree otree;
		/** Reference to Atoms. The Atoms are filled with nulary terms
		 * during execution of parse(). */
		Atoms& atoms;
		/** Vector of formulas (zero derivatives) in the order as they
		 * have been parsed. */
		vector<int> formulas;
		/** The vector to derivatives, each vector corresponds to a
		 * formula in the vector formulas. */
		vector<FormulaDerivatives*> ders;
	public:
		/** Construct an empty formula parser. */
		FormulaParser(Atoms& a)
			: atoms(a) {}
		/** Copy constructor using a different instance of Atoms. */
		FormulaParser(const FormulaParser& fp, Atoms& a);
		virtual ~FormulaParser();

		/** Requires an addition of the formula; called from the
		 * parser. */
		void add_formula(int t);
		/** Requires an addition of the binary operation; called from
		 * the parser. */
		int add_binary(code_t code, int t1, int t2);
		/** Requires an addition of the unary operation; called from
		 * the parser. */
		int add_unary(code_t code, int t);
		/** Requires an addition of the nulary operation given by the
		 * string. The Atoms are consulted for uniquness and are given
		 * an internal index generated by the OperationTree. This is
		 * the channel through which the Atoms are filled. */
		int add_nulary(const char* str);

		/** Adds a derivative to the tree. This just calls
		 * OperationTree::add_derivative. */
		int add_derivative(int t, int v)
			{return otree.add_derivative(t, v);}
		/** Adds a substitution. This just calls
		 * OperationTree::add_substitution. */
		int add_substitution(int t, const map<int,int>& subst)
			{return otree.add_substitution(t, subst);}
		/** Add the substitution given by the map where left sides of
		 * substitutions come from another parser. The right sides are
		 * from this object. The given t is from the given parser fp. */
		int add_substitution(int t, const map<int,int>& subst,
							 const FormulaParser& fp)
			{return otree.add_substitution(t, subst, fp.otree);}
		/** This adds formulas from the given parser with (possibly)
		 * different atoms applying substitutions from the given map
		 * mapping atoms from fp to atoms of the object. */
		void add_subst_formulas(const map<int,int>& subst, const FormulaParser& fp);
		/** Substitute formulas. For each i from 1 through all
		 * formulas, it adds a substitution of the i-th formula and
		 * make it to be i-th formula.*/
		void substitute_formulas(const std::map<int,int>& subst);
		/** This method turns the given term to nulary operation. It
		 * should be used with caution, since this method does not
		 * anything do with atoms, but usually some action is also
		 * needed (at leat to assign the tree index t to some
		 * atom). */
		void nularify(int t)
			{otree.nularify(t);}
		/** Returns a set of nulary terms of the given term. Just
		 * calls OperationTree::nulary_of_term. */
		const hash_set<int>& nulary_of_term(int t) const
			{return otree.nulary_of_term(t);}

		/** Parse a given string containing one or more formulas. The
		 * formulas are parsed and added to the OperationTree and to
		 * the formulas vector. */
		void parse(int length, const char* stream);
		/** Processes a syntax error from bison. */
		void error(const char* mes) const;
		/** Differentiate all the formulas up to the given order. The
		 * variables with respect to which the derivatives are taken
		 * are obtained by Atoms::variables(). If the derivates exist,
		 * they are destroyed and created again (with possibly
		 * different order). */
		void differentiate(int max_order);
		/** Return i-th formula derivatives. */
		const FormulaDerivatives& derivatives(int i) const;

		/** This returns a maximum index of zero derivative formulas
		 * including all nulary terms. This is a mimumum length of the
		 * tree for which it is safe to evaluate zero derivatives of
		 * the formulas. */
		int last_formula() const;
		/** This returns a tree index of the i-th formula in the
		 * vector. */
		int formula(int i) const
			{return formulas[i];}


		/** This returns a tree index of the last formula and pops its
		 * item from the formulas vector. The number of formulas is
		 * then less by one. Returns -1 if there is no formula. If
		 * there are derivatives of the last formula, they are
		 * destroyed and the vector ders is popped from the back. */
		int pop_last_formula();

		/** This returns a number of formulas. */
		int nformulas() const
			{return (int)(formulas.size());}

		/** This returns a reference to atoms. */
		const Atoms& getAtoms() const
			{return atoms;}
		Atoms& getAtoms()
			{return atoms;}
		/** This returns the tree. */
		const OperationTree& getTree() const
			{return otree;}
		OperationTree& getTree()
			{return otree;}

		/** Debug print. */
		void print() const;
	private:
		/** Hide this copy constructor declaration by declaring it as
		 * private. */
		FormulaParser(const FormulaParser& fp);
		/** Destroy all derivatives. */
		void destroy_derivatives();
	};

	/** This is a pure virtual class defining an interface for all
	 * classes which will load the results of formula (zero
	 * derivative) evaluations. A primitive implementation of this
	 * class can be a vector of doubles. */
	class FormulaEvalLoader {
	public:
		virtual ~FormulaEvalLoader() {}
		/** Set the value res for the given formula. The formula is
		 * identified by an index corresponding to the ordering in
		 * which the formulas have been parsed (starting from
		 * zero). */
		virtual void load(int i, double res) = 0;
	};

	/** This class evaluates a selected subset of terms of the
	 * tree. In the protected constructor, one can constraint the
	 * initialization of the evaluation tree to a given number of
	 * terms in the beginning. Using this constructor, one has to make
	 * sure, that the terms in the beginning do not refer to terms
	 * behind the initial part. */
	class FormulaCustomEvaluator {
	protected:
		/** The evaluation tree. */
		EvalTree etree;
		/** The custom tree indices to be evaluated. */
		vector<int> terms;
	public:
		/** Construct from FormulaParser and given list of terms. */
		FormulaCustomEvaluator(const FormulaParser& fp, const vector<int>& ts)
			: etree(fp.otree), terms(ts)
			{}
		/** Construct from OperationTree and given list of terms. */
		FormulaCustomEvaluator(const OperationTree& ot, const vector<int>& ts)
			: etree(ot), terms(ts)
			{}
		/** Evaluate the terms using the given AtomValues and load the
		 * results using the given loader. The loader is called for
		 * each term in the order of the terms. */
		void eval(const AtomValues& av, FormulaEvalLoader& loader);
	protected:
		FormulaCustomEvaluator(const FormulaParser& fp)
			: etree(fp.otree, fp.last_formula()), terms(fp.formulas)
			{}
	};

	/** This class evaluates zero derivatives of the FormulaParser. */
	class FormulaEvaluator : public FormulaCustomEvaluator {
	public:
		/** Construct from FormulaParser. */
		FormulaEvaluator(const FormulaParser& fp)
			: FormulaCustomEvaluator(fp) {}
	};

	/** This is a pure virtual class defining an interface for all
	 * classes which will load the results of formula derivative
	 * evaluations. */
	class FormulaDerEvalLoader {
	public:
		virtual ~FormulaDerEvalLoader() {}
		/** This loads the result of the derivative of the given
		 * order. The semantics of i is the same as in
		 * FormulaEvalLoader::load. The indices of variables with
		 * respect to which the derivative was taken are stored in
		 * memory pointed by vars. These are the tree indices of the
		 * variables. */
		virtual void load(int i, int order, const int* vars, double res) = 0;
	};

	/** This class is a utility class representing the tensor
	 * multindex. It can basically increment itself, and calculate
	 * its offset in the folded tensor. */
	class FoldMultiIndex {
		/** Number of variables. */
		int nvar;
		/** Dimension. */
		int ord;
		/** The multiindex. */
		int* data;
	public:
		/** Initializes to the zero derivative. Order is 0, data is
		 * empty. */
		FoldMultiIndex(int nv);
		/** Initializes the multiindex to zeros or given i. */
		FoldMultiIndex(int nv, int order, int i = 0);
		/** Makes a new multiindex of the same order applying a given
		 * mapping to the indices. The mapping is supposed to be monotone. */
		FoldMultiIndex(int nv, const FoldMultiIndex& mi, const vector<int>& mp);
		/** Shifting constructor. This adds a given number of orders
		 * to the end, copying the last item to the newly added items,
		 * keeping the index ordered. If the index was empty (zero-th
		 * dimension), then zeros are added. */
		FoldMultiIndex(const FoldMultiIndex& fmi, int new_orders);
		/** Copy constructor. */
		FoldMultiIndex(const FoldMultiIndex& fmi);
		/** Desctructor. */
		virtual ~FoldMultiIndex()
			{delete [] data;}
		/** Assignment operator. */
		const FoldMultiIndex& operator=(const FoldMultiIndex& fmi);
		/** Operator < implementing lexicographic ordering within one
		 * order, increasing order across orders. */
		bool operator<(const FoldMultiIndex& fmi) const;
		bool operator==(const FoldMultiIndex& fmi) const;
		/** Increment the multiindex. */
		void increment();
		/** Return offset of the multiindex in the folded tensor. */ 
		int offset() const;
		const int& operator[](int i) const
			{return data[i];}
		/** Return order of the multiindex, i.e. dimension of the
		 * tensor. */ 
		int order() const
			{return ord;}
		/** Return the number of variables. */
		int nv() const
			{return nvar;}
		/** Return the data. */
		const int* ind() const
			{return data;}
		/** Return true if the end of the tensor is reached. The
		 * result of a subsequent increment should be considered
		 * unpredictable. */
		bool past_the_end() const
			{return (ord == 0) || (data[0] == nvar);}
		/** Prints the multiindex in the brackets. */
		void print() const;
	private:
		static int offset_recurse(int* data, int len, int nv);
	};

	/** This class evaluates derivatives of the FormulaParser. */
	class FormulaDerEvaluator {
		/** Its own instance of EvalTree. */
		EvalTree etree;
		/** The indices of derivatives for each formula. This is a
		 * const copy FormulaParser::ders. We do not allocate nor
		 * deallocate anything here. */
		vector<const FormulaDerivatives*> ders;
		/** A copy of tree indices corresponding to atoms to with
		 * respect the derivatives were taken. */
		vector<int> der_atoms;
	public:
		/** Construct the object from FormulaParser. */
		FormulaDerEvaluator(const FormulaParser& fp);
		/** Evaluate the derivatives from the FormulaParser wrt to all
		 * atoms in variables vector at the given AtomValues. The
		 * given loader is used for output. */
		void eval(const AtomValues& av, FormulaDerEvalLoader& loader, int order);
		/** Evaluate the derivatives from the FormulaParser wrt to a
		 * selection of atoms of the atoms in der_atoms vector at the
		 * given AtomValues. The selection is given by a monotone
		 * mapping to the indices (not values) of the der_atoms. */
		void eval(const vector<int>& mp, const AtomValues& av, FormulaDerEvalLoader& loader,
				  int order);
	};
};

#endif

// Local Variables:
// mode:C++
// End:
