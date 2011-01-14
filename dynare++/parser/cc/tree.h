// Copyright (C) 2005-2011, Ondra Kamenik

#ifndef OGP_TREE_H
#define OGP_TREE_H

#include <vector>
#include <set>
#include <map>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <cstdio>

namespace ogp {

	using boost::unordered_set;
	using boost::unordered_map;
	using std::vector;
	using std::set;
	using std::map;

	/** Enumerator representing nulary, unary and binary operation
	 * codes. For nulary, 'none' is used. When one is adding a new
	 * codes, he should update the code of #OperationTree::add_unary,
	 * #OperationTree::add_binary, and of course
	 * #OperationTree::add_derivative. */
	enum code_t {NONE, UMINUS, LOG, EXP, SIN, COS, TAN, SQRT, ERF,
				 ERFC, PLUS, MINUS, TIMES, DIVIDE, POWER};

	/** Class representing a nulary, unary, or binary operation. */
	class Operation {
	protected:
		/** Code of the operation. */
		code_t code;
		/** First operand. If none, then it is -1. */
		int op1;
		/** Second operand. If none, then it is -1. */
		int op2;

	public:
		/** Constructs a binary operation. */
		Operation(code_t cd, int oper1, int oper2)
			: code(cd), op1(oper1), op2(oper2) {}
		/** Constructs a unary operation. */
		Operation(code_t cd, int oper1)
			: code(cd), op1(oper1), op2(-1) {}
		/** Constructs a nulary operation. */
		Operation()
			: code(NONE), op1(-1), op2(-1) {}
		/** A copy constructor. */
		Operation(const Operation& op)
			: code(op.code), op1(op.op1), op2(op.op2) {}

		/** Operator =. */
		const Operation& operator=(const Operation& op)
			{
				code = op.code;
				op1 = op.op1;
				op2 = op.op2;
				return *this;
			}
		/** Operator ==. */
		bool operator==(const Operation& op) const
			{
				return code == op.code && op1 == op.op1 && op2 == op.op2;
			}
		/** Operator < implementing lexicographic ordering. */
		bool operator<(const Operation& op) const
			{
				return (code < op.code ||
						code == op.code &&
						(op1 < op.op1 || op1 == op.op1 && op2 < op.op2));
			}
		/** Returns a number of operands. */
		int nary() const
			{
				return (op2 == -1)? ((op1 == -1) ? 0 : 1) : 2;
			}
		/** Returns a hash value of the operation. */
		size_t hashval() const
			{
				return op2+1 + (op1+1)^15 + code^30;
			}

		code_t getCode() const
			{ return code; }
		int getOp1() const
			{ return op1; }
		int getOp2() const
			{ return op2; }

	};

	/** This struct is a predicate for ordering of the operations in
	 * OperationTree class. now obsolete */
	struct ltoper {
		bool operator()(const Operation& oper1, const Operation& oper2) const
			{return oper1 < oper2;}
	};

	/** Hash function object for Operation. */
	struct ophash {
		size_t operator()(const Operation& op) const
			{ return op.hashval(); }
	};

	/** This struct is a function object selecting some
	 * operations. The operation is given by a tree index. */
	struct opselector {
		virtual bool operator()(int t) const = 0;
		virtual ~opselector() {}
	};

	/** Forward declaration of OperationFormatter. */
	class OperationFormatter;
	class DefaultOperationFormatter;

	/** Forward declaration of EvalTree to make it friend of OperationTree. */
	class EvalTree;

	/** Class representing a set of trees for terms. Each term is
	 * given a unique non-negative integer. The terms are basically
	 * operations whose (integer) operands point to another terms in
	 * the tree. The terms are stored in the vector. Equivalent unary
	 * and binary terms are stored only once. This class guarantees
	 * the uniqueness. The uniqueness of nulary terms is guaranteed by
	 * the caller, since at this level of Operation abstraction, one
	 * cannot discriminate between different nulary operations
	 * (constants, variables). The uniqueness is enforced by the
	 * unordered_map whose keys are operations and values are integers
	 * (indices of the terms).

	 * This class can also make derivatives of a given term with
	 * respect to a given nulary term. I order to be able to quickly
	 * recognize zero derivativates, we maintain a list of nulary
	 * terms contained in the term. A possible zero derivative is then quickly
	 * recognized by looking at the list. The list is implemented as a
	 * unordered_set of integers.
	 *
	 * In addition, many term can be differentiated multiple times wrt
	 * one variable since they can be referenced multiple times. To
	 * avoid this, for each term we maintain a map mapping variables
	 * to the derivatives of the term. As the caller will
	 * differentiate wrt more and more variables, these maps will
	 * become richer and richer.
	 */
	class OperationTree {
		friend class EvalTree;
		friend class DefaultOperationFormatter;
	protected:
		/** This is the vector of the terms. An index to this vector
		 * uniquelly determines the term. */
		vector<Operation> terms;

		/** This defines a type for a map mapping the unary and binary
		 * operations to their indices. */
		typedef unordered_map<Operation, int, ophash> _Topmap;
		typedef _Topmap::value_type _Topval;

		/** This is the map mapping the unary and binary operations to
		 * the indices of the terms.*/
		_Topmap opmap;

		/** This is a type for a set of integers. */
		typedef unordered_set<int> _Tintset;
		/** This is a vector of integer sets corresponding to the
		 * nulary terms contained in the term. */
		vector<_Tintset> nul_incidence;

		/** This is a type of the map from variables (nulary terms) to
		 * the terms. */
		typedef unordered_map<int, int> _Tderivmap;
		/** This is a vector of derivative mappings. For each term, it
		 * maps variables to the derivatives of the term with respect
		 * to the variables. */
		vector<_Tderivmap> derivatives;

		/** The tree index of the last nulary term. */
		int last_nulary;
	public:
		/** This is a number of constants set in the following
		 * enum. This number reserves space in a vector of terms for
		 * the constants. */
		static const int num_constants = 4;
		/** Enumeration for special terms. We need zero, one, nan and
		 * 2/pi.  These will be always first four terms having indices
		 * zero, one and two, three. If adding anything to this
		 * enumeration, make sure you have updated num_constants above.*/
		enum {zero=0, one=1, nan=2, two_over_pi=3};

		/** The unique constructor which initializes the object to
		 * contain only zero, one and nan and two_over_pi.*/
		OperationTree();

		/** Copy constructor. */
		OperationTree(const OperationTree& ot)
			: terms(ot.terms), opmap(ot.opmap), nul_incidence(ot.nul_incidence),
			  derivatives(ot.derivatives),
			  last_nulary(ot.last_nulary)
			{}

		/** Add a nulary operation. The caller is responsible for not
		 * inserting two semantically equivalent nulary operations.
		 * @return newly allocated index
		 */
		int add_nulary();

		/** Add a unary operation. The uniqness is checked, if it
		 * already exists, then it is not added.
		 * @param code the code of the unary operation
		 * @param op the index of the operand
		 * @return the index of the operation
		*/
		int add_unary(code_t code, int op);

		/** Add a binary operation. The uniqueness is checked, if it
		 * already exists, then it is not added. 
		 * @param code the code of the binary operation
		 * @param op1 the index of the first operand
		 * @param op2 the index of the second operand
		 * @return the index of the operation
		 */
		int add_binary(code_t code, int op1, int op2);

		/** Add the derivative of the given term with respect to the
		 * given nulary operation.
		 * @param t the index of the operation being differentiated
		 * @param v the index of the nulary operation
		 * @return the index of the derivative
		 */
		int add_derivative(int t, int v);

		/** Add the substitution given by the map. This adds a new
		 * term which is equal to the given term with applied
		 * substitutions given by the map replacing each term on the
		 * left by a term on the right. We do not check that the terms
		 * on the left are not subterms of the terms on the right. If
		 * so, the substituted terms are not subject of further
		 * substitution. */
		int add_substitution(int t, const map<int,int>& subst);

		/** Add the substitution given by the map where left sides of
		 * substitutions come from another tree. The right sides are
		 * from this tree. The given t is from the given otree. */
		int add_substitution(int t, const map<int,int>& subst,
							 const OperationTree& otree);

		/** This method turns the given term to a nulary
		 * operation. This is an only method, which changes already
		 * existing term (all other methods add something new). User
		 * should use this with caution and must make sure that
		 * something similar has happened for atoms. In addition, it
		 * does not do anything with derivatives, so it should not be
		 * used after some derivatives were created, and derivatives
		 * already created and saved in derivatives mappings should be
		 * forgotten with forget_derivative_maps. */
		void nularify(int t);

		/** Return the set of nulary terms of the given term. */
		const unordered_set<int>& nulary_of_term(int t) const
			{return nul_incidence[t];}

		/** Select subterms of the given term according a given
		 * operation selector and return the set of terms that
		 * correspond to the compounded operations. The given term is
		 * a compound function of the returned subterms and the
		 * function consists only from operations which yield false in
		 * the selector. */
		unordered_set<int> select_terms(int t, const opselector& sel) const;

		/** Select subterms of the given term according a given
		 * operation selector and return the set of terms that
		 * correspond to the compounded operations. The given term is
		 * a compound function of the returned subterms and the
		 * subterms are maximal subterms consisting from operations
		 * yielding true in the selector. */
		unordered_set<int> select_terms_inv(int t, const opselector& sel) const;

		/** This forgets all the derivative mappings. It is used after
		 * a term has been nularified, and then the derivative
		 * mappings carry wrong information. Note that the derivatives
		 * mappings serve only as a tool for quick returns in
		 * add_derivative. Resseting the mappings is harmless, all the
		 * information is rebuilt in add_derivative without any
		 * additional nodes (trees). */
		void forget_derivative_maps();

		/** This returns an operation of a given term. */
		const Operation& operation(int t) const
			{return terms[t];}

		/** This outputs the operation to the given file descriptor
		 * using the given OperationFormatter. */
		void print_operation_tree(int t, FILE* fd, OperationFormatter& f) const;

		/** Debug print of a given operation: */
		void print_operation(int t) const;

		/** Return the last tree index of a nulary term. */
		int get_last_nulary() const
			{return last_nulary;}

		/** Get the number of all operations. */
		int get_num_op() const
			{return (int)(terms.size());}
	private:
		/** This registers a calculated derivative of the term in the
		 * #derivatives vector.
		 * @param t the index of the term for which we register the derivative
		 * @param v the index of the nulary term (variable) to which
		 * respect the derivative was taken
		 * @param tder the index of the resulting derivative
		 */
		void register_derivative(int t, int v, int tder);
		/** This does the same job as select_terms with the only
		 * difference, that it adds the terms to the given set and
		 * hence can be used recursivelly. */
		void select_terms(int t, const opselector& sel, unordered_set<int>& subterms) const; 
		/** This does the same job as select_terms_inv with the only
		 * difference, that it adds the terms to the given set and
		 * hence can be used recursivelly and returns true if the term
		 * was selected. */
		bool select_terms_inv(int t, const opselector& sel, unordered_set<int>& subterms) const; 
		/** This updates nul_incidence information after the term t
		 * was turned to a nulary term in all terms. It goes through
		 * the tree from simplest terms to teh more complex ones and
		 * changes the nul_incidence information where necesary. It
		 * maintains a set where the changes have been made.*/
		void update_nul_incidence_after_nularify(int t);
	};

	/** EvalTree class allows for an evaluation of the given tree for
	 * a given values of nulary terms. For each term in the
	 * OperationTree the class maintains a resulting value and a flag
	 * if the value has been calculated or set. The life cycle of the
	 * class is the following: After it is initialized, the user must
	 * set values for necessary nulary terms. Then the object can be
	 * requested to evaluate particular terms. During this process,
	 * the number of evaluated terms is increasing. Then the user can
	 * request overall reset of evaluation flags, set the nulary terms
	 * to new values and evaluate a number of terms.
	 *
	 * Note that currently the user cannot request a reset of
	 * evaluation flags only for those terms depending on a given
	 * nulary term. This might be added in future and handeled by a
	 * subclasses of OperationTree and EvalTree, since we need a
	 * support for this in OperationTree.
	 */
	class EvalTree {
	protected:
		/** Reference to the OperationTree over which all evaluations
		 * are done. */
		const OperationTree& otree;
		/** The array of values. */
		double* const values;
		/** The array of evaluation flags. */
		bool* const flags;
		/** The index of last operation in the EvalTree. Length of
		 * values and flags will be then last_operation+1. */
		int last_operation;
	public:
		/** Initializes the evaluation tree for the given operation
		 * tree. If last is greater than -1, that the evaluation tree
		 * will contain only formulas up to the given last index
		 * (included). */
		EvalTree(const OperationTree& otree, int last = -1);
		virtual ~EvalTree()
			{ delete [] values; delete [] flags; }
		/** Set evaluation flag to all terms (besides the first
		 * special terms) to false. */
		void reset_all();
		/** Set value for a given nulary term. */
		void set_nulary(int t, double val);
		/** Evaluate the given term with nulary terms set so far. */
		double eval(int t);
		/** Debug print. */
		void print() const;
		/* Return the operation tree. */
		const OperationTree& getOperationTree() const
			{return otree;}
	private:
		EvalTree(const EvalTree&);
	};

	/** This is an interface describing how a given operation is
	 * formatted for output. */
	class OperationFormatter {
	public:
		/** Empty virtual destructor. */
		virtual ~OperationFormatter() {}
		/** Print the formatted operation op with a given tree index t
		 * to a given descriptor. (See class OperationTree to know
		 * what is a tree index.) This prints all the tree. This
		 * always writes equation, left hand side is a string
		 * represenation (a variable, temporary, whatever) of the
		 * term, the right hand side is a string representation of the
		 * operation (which will refer to other string representation
		 * of subterms). */
		virtual void format(const Operation& op, int t, FILE* fd)=0;
	};

	/** The default formatter formats the formulas with a usual syntax
	 * (for example Matlab). A formatting of atoms and terms might be
	 * reimplemented by a subclass. In addition, during its life, the
	 * object maintains a set of tree indices which have been output
	 * and they are not output any more. */
	class DefaultOperationFormatter : public OperationFormatter {
	protected:
		const OperationTree& otree;
		set<int> stop_set;
	public:
		DefaultOperationFormatter(const OperationTree& ot)
			: otree(ot) {}
		/** Format the operation with the default syntax. */
		void format(const Operation& op, int t, FILE* fd);
		/** This prints a string represenation of the given term, for
		 * example 'tmp10' for term 10. In this implementation it
		 * prints $10. */
		virtual void format_term(int t, FILE* fd) const;
		/** Print a string representation of the nulary term. */
		virtual void format_nulary(int t, FILE* fd) const;
		/** Print a delimiter between two statements. By default it is
		 * "\n". */
		virtual void print_delim(FILE* fd) const;
	};

	class NularyStringConvertor {
	public:
		virtual ~NularyStringConvertor() {}
		/** Return the string representation of the atom with the tree
		 * index t. */
		virtual std::string convert(int t) const = 0;
	};

	/** This class converts the given term to its mathematical string representation. */
	class OperationStringConvertor {
	protected:
		const NularyStringConvertor& nulsc;
		const OperationTree& otree;
	public:
		OperationStringConvertor(const NularyStringConvertor& nsc, const OperationTree& ot)
			: nulsc(nsc), otree(ot) {}
		/** Empty virtual destructor. */
		virtual ~OperationStringConvertor() {}
		/** Convert the operation to the string mathematical
		 * representation. This does not write any equation, just
		 * returns a string representation of the formula. */
		std::string convert(const Operation& op, int t) const;
	};
};

#endif

// Local Variables:
// mode:C++
// End:
