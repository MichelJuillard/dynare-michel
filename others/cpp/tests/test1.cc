#include <iostream>
#include "dynare_cpp_driver.hh"
#include "DecisionRules.hh"

DynareInfo *preprocessorOutput(void);
void steadystate(const double *, const double *, double *, int *);
void FirstDerivatives(const double *y, double *x, int nb_row_x, double *params, double *steady_state, int it_, double *residual, double *g1, double *v2, double *v3);

main(int argc, char **argv)
{
  DynareInfo model_info;

  int endo_nbr = model_info.get_endo_nbr();
  int exo_nbr = model_info.get_exo_nbr();
  double *params = model_info.get_params_data();
  // Steady state
  double *steady_state = new double[endo_nbr];
  int info;
  steadystate(NULL,params, steady_state, &info);
  for (int i=0; i < endo_nbr; ++i)
    std::cout << model_info.get_endo_name_by_index(i) << " " << steady_state[i] << "\n";

  // 1st order approximation
  double qz_criterium = 1.000001;

  vector<size_t> zeta_back = model_info.get_zeta_back();
  vector<size_t> zeta_fwrd = model_info.get_zeta_fwrd();
  vector<size_t> zeta_mixed = model_info.get_zeta_mixed();
  vector<size_t> zeta_static = model_info.get_zeta_static();
  int nfwrd = zeta_fwrd.size();
  int nback = zeta_back.size();
  int nmixed = zeta_mixed.size();
  int nstatic = zeta_static.size();

  int jacob_cols = zeta_back.size()+zeta_fwrd.size()+zeta_mixed.size()+exo_nbr;

  double *exo_steady_state = new double[exo_nbr];
  double *jacob_data = new double[endo_nbr*jacob_cols];

  FirstDerivatives(steady_state, exo_steady_state, 0, params, steady_state, 1, NULL, jacob_data, NULL, NULL);

  DecisionRules dr(endo_nbr, exo_nbr, zeta_fwrd, zeta_back, zeta_mixed,
                   zeta_static, qz_criterium);

  int sdyn = nfwrd + nback + 2*nmixed;
  int jacobian_col_nbr = sdyn + endo_nbr + exo_nbr;
  MatrixView jacob_tmp(jacob_data, endo_nbr, jacobian_col_nbr, endo_nbr);

  Matrix jacobian(endo_nbr, jacobian_col_nbr), g_y(endo_nbr, nback+nmixed), g_u(endo_nbr, exo_nbr);
  jacobian = jacob_tmp;

  try
    {
      dr.compute(jacobian, g_y, g_u);
    }
  catch (GeneralizedSchurDecomposition::GSDException &e)
    {
      std::cerr << e << std::endl;
    }
  catch (DecisionRules::BlanchardKahnException &e)
    {
      std::cerr << e << std::endl;
    }

  Vector eig_real(sdyn), eig_cmplx(sdyn);
  dr.getGeneralizedEigenvalues(eig_real, eig_cmplx);
  std::cout << "Eigenvalues (real part): " << eig_real
            << "Eigenvalues (complex part): " << eig_cmplx << std::endl
            << "g_y = " << std::endl << g_y << std::endl
            << "g_u = " << std::endl << g_u;
}
