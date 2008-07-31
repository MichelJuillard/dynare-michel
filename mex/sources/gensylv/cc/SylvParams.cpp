/*
 * Copyright (C) 2003-2008 Ondra Kamenik
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvParams.cpp,v 1.1.1.1 2004/06/04 13:00:52 kamenik Exp $ */

/* Tag $Name:  $ */

#include "SylvParams.h"

#ifdef MWTYPES_NOT_DEFINED
typedef int mwSize;
#endif

void SylvParams::print(const char* prefix) const
{
	print(stdout, prefix);
}

void SylvParams::print(FILE* fdesc, const char* prefix) const
{
	rcondA1.print(fdesc, prefix,  "reci. cond1 A      ", "%8.4g");
	rcondAI.print(fdesc, prefix,  "reci. condInf A    ", "%8.4g");
	bs_norm.print(fdesc, prefix,  "log10 diag norm    ", "%8.4g");
	f_err1.print(fdesc, prefix,   "abs. err 1 F diag  ", "%8.4g");
	f_errI.print(fdesc, prefix,   "abs. err I F diag  ", "%8.4g");
	viv_err1.print(fdesc, prefix, "abs. err 1 V*invV  ", "%8.4g");
	viv_errI.print(fdesc, prefix, "abs. err I V*invV  ", "%8.4g");
	ivv_err1.print(fdesc, prefix, "abs. err 1 invV*V  ", "%8.4g");
	ivv_errI.print(fdesc, prefix, "abs. err I invV*V  ", "%8.4g");
	f_blocks.print(fdesc, prefix, "num blocks in F    ", "%d");
	f_largest.print(fdesc, prefix,"largest block in F ", "%d");
	f_zeros.print(fdesc, prefix,  "num zeros in F     ", "%d");
	f_offdiag.print(fdesc, prefix,"num offdiag in F   ", "%d");
	if (*method == iter) {
		converged.print(fdesc, prefix,       "converged          ", "%d");
		convergence_tol.print(fdesc, prefix, "convergence tol.   ", "%8.4g");
		iter_last_norm.print(fdesc, prefix,  "last norm          ", "%8.4g");
		max_num_iter.print(fdesc, prefix,    "max num iter       ", "%d");
		num_iter.print(fdesc, prefix,        "num iter           ", "%d");
	} else {
		eig_min.print(fdesc, prefix,         "minimum eigenvalue ", "%8.4g");
	}
	mat_err1.print(fdesc, prefix, "rel. matrix norm1  ", "%8.4g");
	mat_errI.print(fdesc, prefix, "rel. matrix normInf", "%8.4g");
	mat_errF.print(fdesc, prefix, "rel. matrix normFro", "%8.4g");
	vec_err1.print(fdesc, prefix, "rel. vector norm1  ", "%8.4g");
	vec_errI.print(fdesc, prefix, "rel. vector normInf", "%8.4g");
	cpu_time.print(fdesc, prefix, "time (CPU secs)    ", "%8.4g");
}

void SylvParams::copy(const SylvParams& p)
{
	method = p.method;
	convergence_tol = p.convergence_tol;
	max_num_iter = p.max_num_iter;
	bs_norm = p.bs_norm;
	want_check = p.want_check;
	converged = p.converged;
	iter_last_norm = p.iter_last_norm;
	num_iter = p.num_iter;
	f_err1 = p.f_err1;
	f_errI = p.f_errI;
	viv_err1 = p.viv_err1;
	viv_errI = p.viv_errI;
	ivv_err1 = p.ivv_err1;
	ivv_errI = p.ivv_errI;
	f_blocks = p.f_blocks;
	f_largest = p.f_largest;
	f_zeros = p.f_zeros;
	f_offdiag = p.f_offdiag;
	rcondA1 = p.rcondA1;
	rcondAI = p.rcondAI;
	eig_min = p.eig_min;
	mat_err1 = p.mat_err1;
	mat_errI = p.mat_errI;
	mat_errF = p.mat_errF;
	vec_err1 = p.vec_err1;
	vec_errI = p.vec_errI;
	cpu_time = p.cpu_time;
}

void SylvParams::setArrayNames(int& num, const char** names) const
{
	num = 0;
	if (method.getStatus() != undef)
		names[num++] = "method";
	if (convergence_tol.getStatus() != undef)
		names[num++] = "convergence_tol";
	if (max_num_iter.getStatus() != undef)
		names[num++] = "max_num_iter";
	if (bs_norm.getStatus() != undef)
		names[num++] = "bs_norm";
	if (converged.getStatus() != undef)
		names[num++] = "converged";
	if (iter_last_norm.getStatus() != undef)
		names[num++] = "iter_last_norm";
	if (num_iter.getStatus() != undef)
		names[num++] = "num_iter";
	if (f_err1.getStatus() != undef)
		names[num++] = "f_err1";
	if (f_errI.getStatus() != undef)
		names[num++] = "f_errI";
	if (viv_err1.getStatus() != undef)
		names[num++] = "viv_err1";
	if (viv_errI.getStatus() != undef)
		names[num++] = "viv_errI";
	if (ivv_err1.getStatus() != undef)
		names[num++] = "ivv_err1";
	if (ivv_errI.getStatus() != undef)
		names[num++] = "ivv_errI";
	if (f_blocks.getStatus() != undef)
		names[num++] = "f_blocks";
	if (f_largest.getStatus() != undef)
		names[num++] = "f_largest";
	if (f_zeros.getStatus() != undef)
		names[num++] = "f_zeros";
	if (f_offdiag.getStatus() != undef)
		names[num++] = "f_offdiag";
	if (rcondA1.getStatus() != undef)
		names[num++] = "rcondA1";
	if (rcondAI.getStatus() != undef)
		names[num++] = "rcondAI";
	if (eig_min.getStatus() != undef)
		names[num++] = "eig_min";
	if (mat_err1.getStatus() != undef)
		names[num++] = "mat_err1";
	if (mat_errI.getStatus() != undef)
		names[num++] = "mat_errI";
	if (mat_errF.getStatus() != undef)
		names[num++] = "mat_errF";
	if (vec_err1.getStatus() != undef)
		names[num++] = "vec_err1";
	if (vec_errI.getStatus() != undef)
		names[num++] = "vec_errI";
	if (cpu_time.getStatus() != undef)
		names[num++] = "cpu_time";
}

#ifdef MATLAB
mxArray* SylvParams::DoubleParamItem::createMatlabArray() const
{
    return mxCreateScalarDouble(value);
}

mxArray* SylvParams::IntParamItem::createMatlabArray() const
{
	mxArray* res = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
	*((int*)mxGetData(res)) = value;
	return res;
}

mxArray* SylvParams::BoolParamItem::createMatlabArray() const
{
	if (value)
		return mxCreateString("true");
	else
		return mxCreateString("false");
}

mxArray* SylvParams::MethodParamItem::createMatlabArray() const
{
	if (value == iter)
		return mxCreateString("iterative");
	else
		return mxCreateString("recursive");
}

mxArray* SylvParams::createStructArray() const
{
	const char* names[50];
	int num;
	setArrayNames(num, names);
	const mwSize dims[] = {1, 1};
	mxArray* const res = mxCreateStructArray(2, dims, num, names);

	int i = 0;
	if (method.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, method.createMatlabArray());
	if (convergence_tol.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, convergence_tol.createMatlabArray());
	if (max_num_iter.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, max_num_iter.createMatlabArray());
	if (bs_norm.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, bs_norm.createMatlabArray());
	if (converged.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, converged.createMatlabArray());
	if (iter_last_norm.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, iter_last_norm.createMatlabArray());
	if (num_iter.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, num_iter.createMatlabArray());
	if (f_err1.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, f_err1.createMatlabArray());
	if (f_errI.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, f_errI.createMatlabArray());
	if (viv_err1.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, viv_err1.createMatlabArray());
	if (viv_errI.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, viv_errI.createMatlabArray());
	if (ivv_err1.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, ivv_err1.createMatlabArray());
	if (ivv_errI.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, ivv_errI.createMatlabArray());
	if (f_blocks.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, f_blocks.createMatlabArray());
	if (f_largest.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, f_largest.createMatlabArray());
	if (f_zeros.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, f_zeros.createMatlabArray());
	if (f_offdiag.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, f_offdiag.createMatlabArray());
	if (rcondA1.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, rcondA1.createMatlabArray());
	if (rcondAI.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, rcondAI.createMatlabArray());
	if (eig_min.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, eig_min.createMatlabArray());
	if (mat_err1.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, mat_err1.createMatlabArray());
	if (mat_errI.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, mat_errI.createMatlabArray());
	if (mat_errF.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, mat_errF.createMatlabArray());
	if (vec_err1.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, vec_err1.createMatlabArray());
	if (vec_errI.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, vec_errI.createMatlabArray());
	if (cpu_time.getStatus() != undef)
		mxSetFieldByNumber(res, 0, i++, cpu_time.createMatlabArray());

	return res;
}
#endif

// Local Variables:
// mode:C++
// End:
