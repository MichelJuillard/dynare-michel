/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvParams.h,v 1.1.1.1 2004/06/04 13:00:54 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef SYLV_PARAMS_H
#define SYLV_PARAMS_H

#include <cstdio>
#include <cstring>

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
# include <dynmex.h>
#endif

typedef enum {def, changed, undef} status;

template <class _Type>
struct ParamItem {
protected:
	typedef ParamItem<_Type> _Self;
	status s;
	_Type value;
public:
	ParamItem()
		{s = undef;}
	ParamItem(_Type val)
		{value = val; s = def;}
	ParamItem(const _Self& item)
		{value = item.value; s = item.s;}
	const _Self& operator=(const _Self& item)
		{value = item.value; s = item.s; return *this;}
	const _Self& operator=(const _Type& val)
		{value = val; s = changed; return *this;}
	_Type operator*() const
		{return value;}
	status getStatus() const
		{return s;}
	void print(FILE* f, const char* prefix, const char* str, const char* fmt) const
		{
			if (s == undef)
				return;
			char out[1000];
			strcpy(out, prefix);
			strcat(out, str);
			strcat(out, "= ");
			strcat(out, fmt);
			if (s == def)
				strcat(out, " <default>");
			strcat(out,"\n");
			fprintf(f, out, value);
		} 
};

class SylvParams {
public:
	typedef enum {iter, recurse} solve_method;

protected:
	class DoubleParamItem : public ParamItem<double> {
	public:
		DoubleParamItem() : ParamItem<double>() {}
		DoubleParamItem(double val) : ParamItem<double>(val) {}
		DoubleParamItem(const DoubleParamItem& item) : ParamItem<double>(item) {}
		const DoubleParamItem& operator=(const double& val)
			{ParamItem<double>::operator=(val); return *this;}
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
		mxArray* createMatlabArray() const;
#endif
	};

	class IntParamItem : public ParamItem<int> {
	public:
		IntParamItem() : ParamItem<int>() {}
		IntParamItem(int val) : ParamItem<int>(val) {}
		IntParamItem(const IntParamItem& item) : ParamItem<int>(item) {}
		const IntParamItem& operator=(const int& val)
			{ParamItem<int>::operator=(val); return *this;}
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
		mxArray* createMatlabArray() const;
#endif
	};

	class BoolParamItem : public ParamItem<bool> {
	public:
		BoolParamItem() : ParamItem<bool>() {}
		BoolParamItem(bool val) : ParamItem<bool>(val) {}
		BoolParamItem(const BoolParamItem& item) : ParamItem<bool>(item) {}
		const BoolParamItem& operator=(const bool& val)
			{ParamItem<bool>::operator=(val); return *this;}
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
		mxArray* createMatlabArray() const;
#endif
	};

	class MethodParamItem : public ParamItem<solve_method> {
	public:
		MethodParamItem() : ParamItem<solve_method>() {}
		MethodParamItem(solve_method val) : ParamItem<solve_method>(val) {}
		MethodParamItem(const MethodParamItem& item) : ParamItem<solve_method>(item) {}
		const MethodParamItem operator=(const solve_method& val)
			{ParamItem<solve_method>::operator=(val); return *this;}
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
		mxArray* createMatlabArray() const;
#endif
	};

public:
	// input parameters
	MethodParamItem method; // method of solution: iter/recurse
	DoubleParamItem convergence_tol; // norm for what we consider converged
	IntParamItem max_num_iter; // max number of iterations
	DoubleParamItem bs_norm; // Bavely Stewart log10 of norm for diagonalization
	BoolParamItem want_check; // true => allocate extra space for checks
	// output parameters
	BoolParamItem converged; // true if converged
	DoubleParamItem iter_last_norm; // norm of the last iteration
	IntParamItem num_iter; // number of iterations
	DoubleParamItem f_err1; // norm 1 of diagonalization abs. error C-V*F*inv(V)
	DoubleParamItem f_errI; // norm Inf of diagonalization abs. error C-V*F*inv(V)
	DoubleParamItem viv_err1; // norm 1 of error I-V*inv(V)
	DoubleParamItem viv_errI; // norm Inf of error I-V*inv(V)
	DoubleParamItem ivv_err1; // norm 1 of error I-inv(V)*V
	DoubleParamItem ivv_errI; // norm Inf of error I-inv(V)*V
	IntParamItem f_blocks; // number of diagonal blocks of F
	IntParamItem f_largest; // size of largest diagonal block in F
	IntParamItem f_zeros; // number of off diagonal zeros in F
	IntParamItem f_offdiag; // number of all off diagonal elements in F
	DoubleParamItem rcondA1; // reciprocal cond 1 number of A
	DoubleParamItem rcondAI; // reciprocal cond Inf number of A
	DoubleParamItem eig_min; // minimum eigenvalue of the solved system
	DoubleParamItem mat_err1; // rel. matrix 1 norm of A*X-B*X*kron(C,..,C)-D
	DoubleParamItem mat_errI; // rel. matrix Inf norm of A*X-B*X*kron(C,..,C)-D
	DoubleParamItem mat_errF; // rel. matrix Frob. norm of A*X-B*X*kron(C,..,C)-D
	DoubleParamItem vec_err1; // rel. vector 1 norm of A*X-B*X*kron(C,..,C)-D
	DoubleParamItem vec_errI; // rel. vector Inf norm of A*X-B*X*kron(C,..,C)-D
	DoubleParamItem cpu_time; // time of the job in CPU seconds
	// note: remember to change copy() if adding/removing member

	SylvParams(bool wc = false)
		: method(recurse), convergence_tol(1.e-30), max_num_iter(15),
		  bs_norm(1.3), want_check(wc) {}
	SylvParams(const SylvParams& p)
		{copy(p);}
	const SylvParams& operator=(const SylvParams& p)
		{copy(p); return *this;}
	~SylvParams() {}
	void print(const char* prefix) const;
	void print(FILE* fdesc, const char* prefix) const;
	void setArrayNames(int& num, const char** names) const;
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
	mxArray* createStructArray() const;
#endif
private:
	void copy(const SylvParams& p);
};

#endif /* SYLV_PARAMS_H */


// Local Variables:
// mode:C++
// End:
