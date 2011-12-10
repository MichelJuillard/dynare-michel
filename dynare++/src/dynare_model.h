// Copyright (C) 2005-2011, Ondra Kamenik

#ifndef OGDYN_DYNARE_MODEL
#define OGDYN_DYNARE_MODEL

#include "parser/cc/matrix_parser.h"
#include "parser/cc/atom_assignings.h"

#include "dynare_atoms.h"
#include "twod_matrix.h"

#include "Vector.h"
#include "GeneralMatrix.h"

#include <map>
#include <boost/unordered_set.hpp>

namespace ogdyn {
	using boost::unordered_set;
	using std::map;

	/** This represents an interval in a string by the pair of
	 * positions (including the first, excluding the second). A
	 * position is given by the line and the column within the line
	 * (both starting from 1). */
	struct PosInterval {
		int fl;
		int fc;
		int ll;
		int lc;
		PosInterval() {}
		PosInterval(int ifl, int ifc, int ill, int ilc)
			: fl(ifl), fc(ifc), ll(ill), lc(ilc) {}
		const PosInterval& operator=(const PosInterval& pi)
			{fl = pi.fl; fc = pi.fc; ll = pi.ll; lc = pi.lc; return *this;}
		/** This returns the interval beginning and interval length
		 * within the given string. */
		void translate(const char* beg, int len, const char*& ibeg, int& ilen) const;
		/** Debug print. */
		void print() const
			{printf("fl=%d fc=%d ll=%d lc=%d\n",fl,fc,ll,lc);}
	};

	/** This class is basically a GeneralMatrix but is created from
	 * parsed matrix data. */
	class ParsedMatrix : public TwoDMatrix {
	public:
		/** Construct the object from the parsed data of ogp::MatrixParser. */
		ParsedMatrix(const ogp::MatrixParser& mp);
	};


	class PlannerBuilder;
	class PlannerInfo;
	class ForwSubstBuilder;
	class ForwSubstInfo;
	class MultInitSS;
	class ModelSSWriter;

	/** A subclass is responsible for creating param_vals, init_vals,
	 * and vcov_mat. */
	class DynareModel {
		friend class PlannerBuilder;
		friend class ForwSubstBuilder;
		friend class MultInitSS;
		friend class ModelSSWriter;
	protected:
		/** All atoms for whole model. */
		DynareDynamicAtoms atoms;
		/** Parsed model equations. */
		ogp::FormulaParser eqs;
		/** Order of approximation. */
		int order;
		/** A vector of parameters values created by a subclass. It
		 * is stored with natural ordering (outer) of the parameters
		 * given by atoms. */
		Vector* param_vals;
		/** A vector of initial values created by a subclass. It is
		 * stored with internal ordering given by atoms. */
		Vector* init_vals;
		/** A matrix for vcov. It is created by a subclass. */
		TwoDMatrix* vcov_mat;
		/** Tree index of the planner objective. If there was no
		 * planner objective keyword, the value is set to -1. */
		int t_plobjective;
		/** Tree index of the planner discount. If there was no
		 * planner discount keyword, the value is set to -1. */
		int t_pldiscount;
		/** Pointer to PlannerBuilder, which is created only if the
		 * planner's FOC are added to the model. */
		PlannerBuilder* pbuilder;
		/** Pointer to an object which builds auxiliary variables and
		 * equations to rewrite a model containing multiple leads to
		 * an equivalent model having only +1 leads. */
		ForwSubstBuilder* fbuilder;
		/** Pointer to AtomSubstitutions which are created when the
		 * atoms are being substituted because of multiple lags
		 * etc. It uses also an old copy of atoms, which is
		 * created. */
		ogp::AtomSubstitutions* atom_substs;
		/** Pointer to a copy of original atoms before substitutions
		 * took place. */
		ogp::SAtoms* old_atoms;
	public:
		/** Initializes the object to an empty state. */
		DynareModel();
		/** Construct a new deep copy. */
		DynareModel(const DynareModel& dm);
		virtual ~DynareModel();		
		virtual DynareModel* clone() const = 0;
		const DynareDynamicAtoms& getAtoms() const
			{return atoms;}
		const ogp::FormulaParser& getParser() const
			{return eqs;}
		int getOrder() const
			{return order;}
		/** Return the vector of parameter values. */
		const Vector& getParams() const
			{return *param_vals;}
		Vector& getParams()
			{return *param_vals;}
		/** Return the vector of initial values of endo variables. */
		const Vector& getInit() const
			{return *init_vals;} 
		Vector& getInit()
			{return *init_vals;}
		/** Return the vcov matrix. */
		const TwoDMatrix& getVcov() const
			{return *vcov_mat;}
		TwoDMatrix& getVcov()
			{return *vcov_mat;}
		/** Return planner info. */
		const PlannerInfo* get_planner_info() const;
		/** Return forward substitutions info. */
		const ForwSubstInfo* get_forw_subst_info() const;
		/** Return substitutions info. */
		const ogp::SubstInfo* get_subst_info() const;
		/** This sets initial values given in outer ordering. */
		void setInitOuter(const Vector& x);
		/** This returns true if the given term is a function of
		 * hardwired constants, numerical constants and parameters. */
		bool is_constant_term(int t) const;
		/** Debug print. */
		void print() const;
		/** Dump the model to the output stream. This includes
		 * variable declarations, parameter values, model code,
		 * initval, vcov and order. */
		void dump_model(std::ostream& os) const;
	protected:
		/** Adds a name of endogenous, exogenous or a parameter. The
		 * sort is governed by the flag. See dynglob.y for values of
		 * the flag. This is used by a subclass when declaring the
		 * names. */
		void add_name(const char* name, int flag);
		/** This checks the model consistency. Thus includes: number
		 * of endo variables and number of equations, min and max lag
		 * of endogenous variables and occurrrences of exogenous
		 * variables. It throws an exception, if there is a problem. */
		void check_model() const;
		/** This shifts the given variable identified by the tree
		 * index in time. So if the given tree index represents a(+3)
		 * and the tshift is -4, the method returns tree index of the
		 * a(-1). If a(-1) doesn't exist, it is added to the tree. If
		 * it exists, its tree index is returned. If the tree index
		 * doesn't correspond to an endogenous nor exogenous variable,
		 * an exception is thrown. */
		int variable_shift(int t, int tshift);
		/** For the given set of atoms identified by tree indices and
		 * given time shift, this method returns a map mapping each
		 * variable in the given set to its time shifted variable. The
		 * map is passed through the reference and is cleared in the
		 * beginning. */
		void variable_shift_map(const unordered_set<int>& a_set, int tshift,
								map<int,int>& s_map);
		/** This returns maximum lead and minimum lag of an endogenous
		 * or exogenous variable in the given term. If there are no
		 * endo or exo variables, than it returns the least integer as
		 * max lead and the greatest integer as min lag. */ 
		void termspan(int t, int& mlead, int& mlag) const;
		/** This function returns a set of non-linear subterms of the
		 * given term, these are terms whose linear combination
		 * constitutes the given term. */
		unordered_set<int> get_nonlinear_subterms(int t) const;
		/** This method assigns already used tree index of some term
		 * to the not-yet used atom name with the given lead/lag. In
		 * this way, all occurrences of term t are substituted with
		 * the atom name(ll). The method handles also rewriting
		 * operation tree including derivatives of the term t. */
		void substitute_atom_for_term(const char* name, int ll, int t);
		/** This performs a final job after the model is parsed. It
		 * creates the PlannerBuilder object if the planner's FOC are
		 * needed, then it creates ForwSubstBuilder handling multiple
		 * leads and finally it creates the substitution object saving
		 * old atoms and performs the substitutions. */
		void final_job();
	};

	/** This class constructs DynareModel from dynare++ model file. It
	 * parses variable declarations, model equations, parameter
	 * assignments, initval assignments, vcov matrix and order of
	 * approximation. */
	class DynareParser : public DynareModel {
	protected:
		/** Static atoms for parameter assignments. */
		DynareStaticAtoms pa_atoms;
		/** Assignments for the parameters. */
		ogp::AtomAssignings paramset;
		/** Static atoms for initval assignments. */
		DynareStaticAtoms ia_atoms;
		/** Assignments for the initval. */
		ogp::AtomAssignings initval;
		/** Matrix parser for vcov. */
		ogp::MatrixParser vcov;
	public:
		/** This, in fact, creates DynareModel from the given string
		 * of the given length corresponding to the Dynare++ model
		 * file. If the given ord is not -1, then it overrides setting
		 * in the model file. */
		DynareParser(const char* str, int len, int ord);
		DynareParser(const DynareParser& p);
		virtual ~DynareParser();
		DynareModel* clone() const
			{return new DynareParser(*this);}
		/** Adds a name of endogenous, exogenous or a parameter. This
		 * addss the name to the parent class DynareModel and also
		 * registers the name to either paramset, or initval. */
		void add_name(const char* name, int flag);
		/** Sets position of the model section. Called from
		 * dynglob.y. */
		void set_model_pos(int off1, int off2)
			{model_beg = off1; model_end = off2;}
		/** Sets position of the section setting parameters. Called
		 * from dynglob.y. */
		void set_paramset_pos(int off1, int off2)
			{paramset_beg = off1; paramset_end = off2;}
		/** Sets position of the initval section. Called from
		 * dynglob.y. */
		void set_initval_pos(int off1, int off2)
			{initval_beg = off1; initval_end = off2;}
		/** Sets position of the vcov section. Called from
		 * dynglob.y. */
		void set_vcov_pos(int off1, int off2)
			{vcov_beg = off1; vcov_end = off2;}
		/** Parser the given string as integer and set to as the
		 * order. */
		void set_order_pos(int off1, int off2)
			{order_beg = off1; order_end = off2;}
		/** Sets position of the planner_objective section. Called
		 * from dynglob.y. */
		void set_pl_objective_pos(int off1, int off2)
			{plobjective_beg = off1; plobjective_end = off2;}
		/** Sets position of the planner_discount section. Called from
		 * dynglob.y. */
		void set_pl_discount_pos(int off1, int off2)
			{pldiscount_beg = off1; pldiscount_end = off2;}		
		/** Processes a syntax error from bison. */
		void error(const char* mes);
		/** Debug print. */
		void print() const;
	protected:
		void parse_glob(int length, const char* stream);
		int parse_order(int length, const char* stream);
		int parse_pldiscount(int length, const char* stream);
		/** Evaluate paramset assignings and set param_vals. */
		void calc_params();
		/** Evaluate initval assignings and set init_vals. */
		void calc_init();
		/** Do the final job. This includes building the planner
		 * problem (if any) and substituting for multiple lags, and
		 * one period leads of exogenous variables, and calculating
		 * initial guess of lagrange multipliers in the social planner
		 * problem. Precondtion: everything parsed and calculated
		 * parameters, postcondition: calculated initvals vector and
		 * parsing_finished for expanded vectors. */
		void final_job();
	private:
		int model_beg, model_end;
		int paramset_beg, paramset_end;
		int initval_beg, initval_end;
		int vcov_beg, vcov_end;
		int order_beg, order_end;
		int plobjective_beg, plobjective_end;
		int pldiscount_beg, pldiscount_end;
	};

	/** Semiparsed model. The equations are given by a string,
	 * everything other by C/C++ objects. The initial values are set
	 * manually after the creation of this object. This implies that
	 * no automatic substitutions cannot be done here, which in turn
	 * implies that we cannot do here a social planner nor substitutions
	 * of multiple lags. */
	class DynareSPModel : public DynareModel {
	public:
		DynareSPModel(const char** endo, int num_endo,
					  const char** exo, int num_exo,
					  const char** par, int num_par,
					  const char* equations, int len, int ord);
		DynareSPModel(const DynareSPModel& dm)
			: DynareModel(dm) {}
		~DynareSPModel() {}
		virtual DynareModel* clone() const
			{return new DynareSPModel(*this);}
	};

	/** This class implements a selector of operations which correspond
	 * to non-linear functions. This inherits from ogp::opselector and
	 * is used to calculate non-linear subterms in
	 * DynareModel::get_nonlinear_subterms(). */
	class NLSelector : public ogp::opselector {
	private:
		const DynareModel& model;
	public:
		NLSelector(const DynareModel& m) : model(m) {}
		bool operator()(int t) const;
	};

	/** This class writes a mathematical code evaluating the system of
	 * equations and the first derivatives at zero shocks and at the
	 * given (static) state. Static means that lags and leads are
	 * ignored. */
	class ModelSSWriter : public ogp::DefaultOperationFormatter {
	protected:
		const DynareModel& model;
	public:
		ModelSSWriter(const DynareModel& m)
			: DefaultOperationFormatter(m.eqs.getTree()),
			  model(m) {}
		/** This writes the evaluation of the system. It calls pure
		 * virtual methods for writing a preamble, then assignment of
		 * atoms, and then assignment for resulting object. These are
		 * language dependent and are implemented in the subclass. */
		void write_der0(FILE* fd);
		/** This writes the evaluation of the first order derivative of
		the system. It calls pure virtual methods for writing a
		preamble, assignment, and assignemnt of the resulting
		objects. */
		void write_der1(FILE* fd);
	protected:
		virtual void write_der0_preamble(FILE* fd) const =0;
		virtual void write_der1_preamble(FILE* fd) const =0;
		virtual void write_atom_assignment(FILE* fd) const =0;
		virtual void write_der0_assignment(FILE* fd) const =0;
		virtual void write_der1_assignment(FILE* fd) const =0;
	};


	class MatlabSSWriter : public ModelSSWriter {
	protected:
		/** Identifier used in function names. */
		char* id;
	public:
		MatlabSSWriter(const DynareModel& dm, const char* idd);
		virtual ~MatlabSSWriter()
			{delete [] id;}
	protected:
		// from ModelSSWriter
		void write_der0_preamble(FILE* fd) const;
		void write_der1_preamble(FILE* fd) const;
		/** This writes atom assignments. We have four kinds of atoms
		 * set here: endogenous vars coming from one parameter,
		 * parameter values given by the second parameter, constants,
		 * and the OperationTree::num_constants hardwired constants in
		 * ogp::OperationTree. */
		void write_atom_assignment(FILE* fd) const;
		void write_der0_assignment(FILE* fd) const;
		void write_der1_assignment(FILE* fd) const;
		/** This prints t10 for t=10. */
		void format_term(int t, FILE* fd) const;
		/** This prints a10 for t=10. The atoms a10 are supposed to be
		 * set by write_atom_assignments(). */
		void format_nulary(int t, FILE* fd) const;
	private:
		void write_common1_preamble(FILE* fd) const;
		void write_common2_preamble(FILE* fd) const;
	};

	/** This class implements OperationFormatter for debugging
	 * purposes. It renders atoms in a more friendly way than the
	 * ogp::DefaulOperationFormatter. */
	class DebugOperationFormatter : public ogp::DefaultOperationFormatter {
	protected:
		const DynareModel& model;
	public:
		DebugOperationFormatter(const DynareModel& m)
			: DefaultOperationFormatter(m.getParser().getTree()),
			  model(m) {}
		void format_nulary(int t, FILE* fd) const;
	};
};

#endif

// Local Variables:
// mode:C++
// End:
