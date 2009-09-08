// Copyright (C) 2006, Ondra Kamenik

// $Id$

#include "planner_builder.h"
#include "dynare_exception.h"

#include <cmath>

using namespace ogdyn;

const IntegerMatrix& IntegerMatrix::operator=(const IntegerMatrix& im)
{
	if (nr != im.nr || nc != im.nc)
		throw DynareException(__FILE__,__LINE__,
							  "Matrices have different dimensions in IntegerMatrix::operator=");
	memcpy(data, im.data, nr*nc*sizeof(int));
	return *this;
}

const IntegerArray3& IntegerArray3::operator=(const IntegerArray3& ia3)
{
	if (n1 != ia3.n1 || n2 != ia3.n2 || n3 != ia3.n3)
		throw DynareException(__FILE__,__LINE__,
							  "Arrays have different dimensions in IntegerArray3::operator=");
	memcpy(data, ia3.data, n1*n2*n3*sizeof(int));
	return *this;
}


PlannerBuilder::PlannerBuilder(DynareModel& m, const Tvarset& yyset,
							   const Teqset& ffset)
	: yset(), fset(ffset), model(m),
	  tb(model.t_plobjective), tbeta(model.t_pldiscount),
	  maxlead(model.atoms.get_maxlead()),
	  minlag(model.atoms.get_minlag()),
	  diff_b(yyset.size(), 1-minlag),
	  diff_f(yyset.size(), fset.size(), 1+maxlead-minlag),
	  static_atoms(),
	  static_tree(),
	  diff_b_static(yyset.size(), 1-minlag),
	  diff_f_static(yyset.size(), fset.size(), 1+maxlead-minlag)
{
	info.num_new_terms -= model.getParser().getTree().get_num_op();

	fill_yset(m.atoms.get_name_storage(), yyset);

	add_derivatives_of_b();
	add_derivatives_of_f();
	shift_derivatives_of_b();
	shift_derivatives_of_f();
	beta_multiply_b();
	beta_multiply_f();
	make_static_version();
	lagrange_mult_f();
	form_equations();

	info.num_new_terms += model.getParser().getTree().get_num_op();
}

PlannerBuilder::PlannerBuilder(const PlannerBuilder& pb, ogdyn::DynareModel& m)
	: yset(), fset(pb.fset), model(m),
	  tb(pb.tb), tbeta(pb.tbeta),
	  maxlead(pb.maxlead), minlag(pb.minlag),
	  diff_b(pb.diff_b), diff_f(pb.diff_f),
	  static_atoms(pb.static_atoms),
	  static_tree(pb.static_tree),
	  diff_b_static(pb.diff_b_static),
	  diff_f_static(pb.diff_f_static),
	  aux_map(), static_aux_map()
	
{
	fill_yset(m.atoms.get_name_storage(), pb.yset);
	fill_aux_map(m.atoms.get_name_storage(), pb.aux_map, pb.static_aux_map);
}

void PlannerBuilder::add_derivatives_of_b()
{
	int yi = 0;
	for (Tvarset::const_iterator yname = yset.begin();
		 yname != yset.end(); ++yname, yi++)
		for (int ll = minlag; ll <= 0; ll++) {
			int yt = model.atoms.index(*yname, ll);
			if (yt != -1)
				diff_b(yi, ll-minlag) = model.eqs.add_derivative(tb, yt);
			else
				diff_b(yi, ll-minlag) = ogp::OperationTree::zero;
		}
}

void PlannerBuilder::add_derivatives_of_f()
{
	int yi = 0;
	for (Tvarset::const_iterator yname = yset.begin();
		 yname != yset.end(); ++yname, yi++)
		for (unsigned int fi = 0; fi < fset.size(); fi++)
			for (int ll = minlag; ll <= maxlead; ll++) {
				int yt = model.atoms.index(*yname, ll);
				if (yt != -1)
					diff_f(yi, fi, ll-minlag) =
						model.eqs.add_derivative(model.eqs.formula(fset[fi]), yt);
				else
					diff_f(yi, fi, ll-minlag) = ogp::OperationTree::zero;
			}
}

void PlannerBuilder::shift_derivatives_of_b()
{
	map<int,int> subst;
	for (int yi = 0; yi < diff_b.nrows(); yi++)
		for (int ll = minlag; ll < 0; ll++)
			if (diff_b(yi, ll-minlag) != ogp::OperationTree::zero) {
				model.variable_shift_map(model.eqs.nulary_of_term(diff_b(yi, ll-minlag)),
										 -ll, subst);
				diff_b(yi, ll-minlag) = model.eqs.add_substitution(diff_b(yi, ll-minlag), subst);
			}
}

void PlannerBuilder::shift_derivatives_of_f()
{
	map<int,int> subst;
	for (int yi = 0; yi < diff_f.dim1(); yi++)
		for (int fi = 0; fi < diff_f.dim2(); fi++) {
			// first do it leads which are put under expectation before t: no problem
			for (int ll = 0; ll <= maxlead; ll++)
				if (diff_f(yi, fi, ll-minlag) != ogp::OperationTree::zero) {
					model.variable_shift_map(model.eqs.nulary_of_term(diff_f(yi, fi, ll-minlag)),
											 -ll, subst);
					diff_f(yi, fi, ll-minlag) =
						model.eqs.add_substitution(diff_f(yi, fi, ll-minlag), subst);
				}
			// now do it for lags, these are put as leads under
			// expectations after time t, so we have to introduce
			// auxiliary variables at time t, and make leads of them here
			for (int ll = minlag; ll < 0; ll++) {
				int ft = diff_f(yi, fi, ll-minlag);
				if (ft != ogp::OperationTree::zero) {
					// if the ft term has a lead, than we need to
					// introduce an auxiliary variable z_t, define it
					// as E_t[ft] and put z_{t-ll} to the
					// equation. Otherwise, we just put leaded ft to
					// the equation directly.
					int ft_maxlead, ft_minlag;
					model.termspan(ft, ft_maxlead, ft_minlag);
					if (ft_maxlead > 0) {
						// make an auxiliary variable
						char name[100];
						sprintf(name, "AUX_%d_%d_%d", yi, fset[fi], -ll);
						model.atoms.register_uniq_endo(name);
						info.num_aux_variables++;
						int taux = model.eqs.add_nulary(name);
						sprintf(name, "AUX_%d_%d_%d(%d)", yi, fset[fi], -ll, -ll);
						int taux_leaded = model.eqs.add_nulary(name);
						// put aux_leaded to the equation
						diff_f(yi, fi, ll-minlag) = taux_leaded;
						// save auxiliary variable and the term
						aux_map.insert(Tsubstmap::value_type(model.atoms.name(taux), ft));
					} else {
						// no auxiliary variable is needed and the
						// term ft can be leaded in place
						model.variable_shift_map(model.eqs.nulary_of_term(ft), -ll, subst);
						diff_f(yi, fi, ll-minlag) =
							model.eqs.add_substitution(ft, subst);
					}
				}
			}
		}
}

void PlannerBuilder::beta_multiply_b()
{
	int beta_pow = ogp::OperationTree::one;
	for (int ll = 0; ll >= minlag; ll--,
			 beta_pow = model.eqs.add_binary(ogp::TIMES, beta_pow, tbeta))
		for (int yi = 0; yi < diff_b.nrows(); yi++)
			if (diff_b(yi, ll-minlag) != ogp::OperationTree::zero)
				diff_b(yi, ll-minlag) =
					model.eqs.add_binary(ogp::TIMES, beta_pow, diff_b(yi, ll-minlag));
}

void PlannerBuilder::beta_multiply_f()
{
	int beta_pow = ogp::OperationTree::one;
	for (int ll = 0; ll <= maxlead; ll++,
			 beta_pow = model.eqs.add_binary(ogp::DIVIDE, beta_pow, tbeta))
		for (int yi = 0; yi < diff_f.dim1(); yi++)
			for (int fi = 0; fi < diff_f.dim2(); fi++)
				if (diff_f(yi, fi, ll-minlag) != ogp::OperationTree::zero)
					diff_f(yi, fi, ll-minlag) =
						model.eqs.add_binary(ogp::TIMES, beta_pow, diff_f(yi, fi, ll-minlag));

	beta_pow = ogp::OperationTree::one;
	for (int ll = 0; ll >= minlag; ll--,
			 beta_pow = model.eqs.add_binary(ogp::TIMES, beta_pow, tbeta))
		for (int yi = 0; yi < diff_f.dim1(); yi++)
			for (int fi = 0; fi < diff_f.dim2(); fi++)
				if (diff_f(yi, fi, ll-minlag) != ogp::OperationTree::zero)
					diff_f(yi, fi, ll-minlag) =
						model.eqs.add_binary(ogp::TIMES, beta_pow, diff_f(yi, fi, ll-minlag));
}

void PlannerBuilder::make_static_version()
{
	// map holding substitutions from dynamic to static
	ogp::StaticFineAtoms::Tintintmap tmap;

	// fill static atoms with outer ordering
	static_atoms.import_atoms(model.atoms, static_tree, tmap);

	// go through diff_b and fill diff_b_static
	for (int ll = minlag; ll <= 0; ll++)
		for (int yi = 0; yi < diff_b.nrows(); yi++)
			diff_b_static(yi, ll-minlag) =
				static_tree.add_substitution(diff_b(yi, ll-minlag),
											 tmap,  model.eqs.getTree());

	// go through diff_f and fill diff_f_static
	for (int ll = minlag; ll <= maxlead; ll++)
		for (int yi = 0; yi < diff_f.dim1(); yi++)
			for (int fi = 0; fi < diff_f.dim2(); fi++)
				diff_f_static(yi, fi, ll-minlag) =
					static_tree.add_substitution(diff_f(yi, fi, ll-minlag),
												 tmap, model.eqs.getTree());

	// go through aux_map and fill static_aux_map
	for (Tsubstmap::const_iterator it = aux_map.begin();
		 it != aux_map.end(); ++it) {
		int tstatic = static_tree.add_substitution((*it).second, tmap, model.eqs.getTree());
		const char* name = static_atoms.get_name_storage().query((*it).first);
		static_aux_map.insert(Tsubstmap::value_type(name, tstatic));
	}
}


void PlannerBuilder::lagrange_mult_f()
{
	// register multipliers
	char mult_name[100];
	for (int fi = 0; fi < diff_f.dim2(); fi++) {	
		sprintf(mult_name, "MULT%d", fset[fi]);
		model.atoms.register_uniq_endo(mult_name);
		info.num_lagrange_mults++;
	}
	// multiply with the multipliers
	for (int yi = 0; yi < diff_f.dim1(); yi++)
		for (int fi = 0; fi < diff_f.dim2(); fi++)
			for (int ll = minlag; ll <= maxlead; ll++)
				if (diff_f(yi, fi, ll-minlag) != ogp::OperationTree::zero) {
					sprintf(mult_name, "MULT%d(%d)", fset[fi], -ll);
					int tm = model.eqs.add_nulary(mult_name);
					diff_f(yi, fi, ll-minlag) =
						model.eqs.add_binary(ogp::TIMES, tm, diff_f(yi, fi, ll-minlag));
				}
}

void PlannerBuilder::form_equations()
{
	// add planner's FOCs
	for (int yi = 0; yi < diff_f.dim1(); yi++) {
		int eq = ogp::OperationTree::zero;
		for (int ll = minlag; ll <= 0; ll++)
			eq = model.eqs.add_binary(ogp::PLUS, eq, diff_b(yi, ll-minlag));
		for (int fi = 0; fi < diff_f.dim2(); fi++)
			for (int ll = minlag; ll <= maxlead; ll++)
				eq = model.eqs.add_binary(ogp::PLUS, eq, diff_f(yi, fi, ll-minlag));
		model.eqs.add_formula(eq);
	}

	// add equations for auxiliary variables
	for (Tsubstmap::const_iterator it = aux_map.begin();
		 it != aux_map.end(); ++it) {
		int t = model.atoms.index((*it).first, 0);
		model.eqs.add_formula(model.eqs.add_binary(ogp::MINUS, t, (*it).second));
	}
}

void PlannerBuilder::fill_yset(const ogp::NameStorage& ns,
							   const PlannerBuilder::Tvarset& yyset)
{
	for (Tvarset::const_iterator it = yyset.begin(); it != yyset.end(); ++it)
		yset.insert(ns.query(*it));
}

void PlannerBuilder::fill_aux_map(const ogp::NameStorage& ns, const Tsubstmap& aaux_map,
								  const Tsubstmap& astatic_aux_map)
{
	// fill aux_map
	for (Tsubstmap::const_iterator it = aaux_map.begin();
		 it != aaux_map.end(); ++it)
		aux_map.insert(Tsubstmap::value_type(ns.query((*it).first), (*it).second));

	// fill static_aux_map
	for (Tsubstmap::const_iterator it = astatic_aux_map.begin();
		 it != astatic_aux_map.end(); ++it)
		static_aux_map.insert(Tsubstmap::value_type(static_atoms.get_name_storage().query((*it).first),
													(*it).second));
}

MultInitSS::MultInitSS(const PlannerBuilder& pb, const Vector& pvals, Vector& yy)
	: builder(pb), b(builder.diff_b_static.nrows()),
	  F(builder.diff_f_static.dim1(), builder.diff_f_static.dim2())
{
	b.zeros();
	F.zeros();

	// first evaluate substitutions (auxiliary variables) from the builder
	ogdyn::DynareStaticSteadySubstitutions dss(builder.model.atoms, builder.static_atoms,
											   builder.static_tree,
											   builder.static_aux_map, pvals, yy);

	// gather all the terms from builder.diff_b_static and
	// builder.diff_f_static to the vector, the ordering is important,
	// since the index of this vector will have to be decoded to the
	// position in b and F.
	vector<int> terms;
	for (int yi = 0; yi < builder.diff_b_static.nrows(); yi++)
		for (int l = 0; l < builder.diff_b_static.ncols(); l++)
			terms.push_back(builder.diff_b_static(yi, l));
	for (int yi = 0; yi < builder.diff_f_static.dim1(); yi++)
		for (int fi = 0; fi < builder.diff_f_static.dim2(); fi++)
			for (int l = 0; l < builder.diff_f_static.dim3(); l++)
				terms.push_back(builder.diff_f_static(yi, fi, l));

	// evaluate the terms, it will call a series of load(i,res), which
	// sum the results through lags/leads to b and F
	DynareStaticSteadyAtomValues dssav(builder.model.atoms, builder.static_atoms, pvals, yy);
	ogp::FormulaCustomEvaluator fe(builder.static_tree, terms);
	fe.eval(dssav, *this);

	// solve overdetermined system b+F*lambda=0 => lambda=-(F^T*F)^{-1}*F^T*b
	GeneralMatrix FtF(F, "transpose", F);
	Vector lambda(builder.diff_f_static.dim2());
	F.multVecTrans(0.0, lambda, -1.0, b);
	ConstGeneralMatrix(FtF).multInvLeft(lambda);

	// take values of lambda and put it to yy
	for (int fi = 0; fi < builder.diff_f_static.dim2(); fi++) {
		char mult_name[100];
		sprintf(mult_name, "MULT%d", builder.fset[fi]);
		int iouter = builder.model.atoms.name2outer_endo(mult_name);
		int iy = builder.model.atoms.outer2y_endo()[iouter];
		if (! std::isfinite(yy[iy]))
			yy[iy] = lambda[fi];

		// go through all substitutions of the multiplier and set them
		// as well
		if (builder.model.atom_substs) {
			const ogp::AtomSubstitutions::Toldnamemap& old2new =
				builder.model.atom_substs->get_old2new();
			const ogp::AtomSubstitutions::Toldnamemap::const_iterator it =
				old2new.find(mult_name);
			if (it != old2new.end()) {
				const ogp::AtomSubstitutions::Tshiftnameset& sset = (*it).second;
				for (ogp::AtomSubstitutions::Tshiftnameset::const_iterator itt = sset.begin();
					 itt != sset.end(); ++itt) {
					const char* newname = (*itt).first;
					int iouter = builder.model.atoms.name2outer_endo(newname);
					int iy = builder.model.atoms.outer2y_endo()[iouter];
					if (! std::isfinite(yy[iy]))
						yy[iy] = lambda[fi];
				}
			}
		}
	}
}

void MultInitSS::load(int i, double res)
{
	// we can afford it, since the evaluator sets res to exact zero if
	// the term is zero
	if (res == 0)
		return;
	// decode i and add to either b or F 
	if (i < builder.diff_b_static.nrows()*builder.diff_b_static.ncols()) {
		// add to b
		b[i / builder.diff_b_static.ncols()] += res;
	} else {
		// add to F
		i -= builder.diff_b_static.nrows()*builder.diff_b_static.ncols();
		int yifi = i / builder.diff_f_static.dim3();
		int yi = yifi / builder.diff_f_static.dim2();
		int fi = yifi % builder.diff_f_static.dim2();
		F.get(yi, fi) += res;
	}
}
