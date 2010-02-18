/*
 * This test of DecisionRules class is based on example1.mod.
 */

/*
 * Copyright (C) 2010 Dynare Team
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

#include "DecisionRules.hh"

int
main(int argc, char **argv)
{
  size_t endo_nbr = 6, exo_nbr = 2;

  std::vector<size_t> zeta_fwrd, zeta_back, zeta_mixed, zeta_static;
  // y and c are purely forward
  zeta_fwrd.push_back(0);
  zeta_fwrd.push_back(1);
  // k and a are purely backward
  zeta_back.push_back(2);
  zeta_back.push_back(3);
  // h is static
  zeta_static.push_back(4);
  // b is both backward and forward
  zeta_mixed.push_back(5);

  double qz_criterium = 1.000001;

  DecisionRules dr(endo_nbr, exo_nbr, zeta_fwrd, zeta_back, zeta_mixed,
                   zeta_static, qz_criterium);

  double jacob_data[] = {
    0.000000000000000,
    0.000000000000000,
    -0.035101010101010,
    -0.975000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    -0.950000000000000,
    -0.025000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    -0.025000000000000,
    -0.950000000000000,
    -0.640000000000000,
    0.000000000000000,
    1.000000000000000,
    -1.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.860681114551094,
    -13.792569659442703,
    0.000000000000000,
    1.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.034750000000000,
    0.000000000000000,
    1.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    -1.080682530956729,
    0.000000000000000,
    1.000000000000000,
    0.000000000000000,
    2.370597639417809,
    0.000000000000000,
    -2.370597639417800,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    -11.083604432603581,
    0.000000000000000,
    -0.277090110815090,
    0.000000000000000,
    1.000000000000000,
    0.000000000000000,
    -0.356400000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    13.792569659442703,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    10.698449178570606,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    -1.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    0.000000000000000,
    -1.000000000000000
  };

  MatrixView jacob_tmp(jacob_data, 6, 14, 6);

  Matrix jacobian(6, 14), g_y(6, 3), g_u(6, 2);
  jacobian = jacob_tmp;

  try
    {
      dr.compute(jacobian, g_y, g_u);
    }
  catch(GeneralizedSchurDecomposition::GSDException &e)
    {
      std::cerr << e << std::endl;
    }
  catch(DecisionRules::BlanchardKahnException &e)
    {
      std::cerr << e << std::endl;
    }

  Vector eig_real(6), eig_cmplx(6);
  dr.getGeneralizedEigenvalues(eig_real, eig_cmplx);
  std::cout << "Eigenvalues (real part): " << eig_real
            << "Eigenvalues (complex part): " << eig_cmplx;
}
