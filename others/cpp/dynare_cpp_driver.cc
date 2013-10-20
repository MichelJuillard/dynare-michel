/*
 * Copyright (C) 2011-2012 Houtan Bastani, Daniel Waggoner, Tao Zha
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * More information is available at <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <assert.h>

using namespace std;

#include "dynare_cpp_driver.hh"

DynareInfo::DynareInfo(map<string, int > exo_names_arg,
                       map<string, int > exo_det_names_arg,
                       map<string, int > endo_names_arg,
                       map<string, int > param_names_arg,
                       vector< double > params_arg,
                       vector<aux_vars_t> aux_vars_arg,
                       vector<int> predetermined_variables_arg,
                       vector<int> varobs_arg,
                       vector<vector<int > > lead_lag_incidence_arg,
                       vector<int> NNZDerivatives_arg) :
  exo_names(exo_names_arg),
  exo_det_names(exo_det_names_arg),
  endo_names(endo_names_arg),
  param_names(param_names_arg),
  params(params_arg),
  aux_vars(aux_vars_arg),
  predetermined_variables(predetermined_variables_arg),
  varobs(varobs_arg),
  lead_lag_incidence(lead_lag_incidence_arg),
  NNZDerivatives(NNZDerivatives_arg)
{
}

DynareInfo::~DynareInfo()
{
  for (vector<MarkovSwitching *>::iterator it = markov_switching_vector.begin();
       it < markov_switching_vector.end(); it++ )
    delete *it;

  for (vector<ModFilePrior *>::iterator it = prior_vector.begin();
       it < prior_vector.end(); it++ )
    delete *it;

  for (vector<ModFileStructuralInnovationPrior *>::iterator it = structural_innovation_prior_vector.begin();
       it < structural_innovation_prior_vector.end(); it++ )
    delete *it;

  for (vector<ModFileMeasurementErrorPrior *>::iterator it = measurement_error_prior_vector.begin();
       it < measurement_error_prior_vector.end(); it++ )
    delete *it;

  for (vector<ModFileStructuralInnovationCorrPrior *>::iterator it = structural_innovation_corr_prior_vector.begin();
       it < structural_innovation_corr_prior_vector.end(); it++ )
    delete *it;

  for (vector<ModFileMeasurementErrorCorrPrior *>::iterator it = measurement_error_corr_prior_vector.begin();
       it < measurement_error_corr_prior_vector.end(); it++ )
    delete *it;

  markov_switching_vector.clear();
  prior_vector.clear();
  structural_innovation_prior_vector.clear();
  measurement_error_prior_vector.clear();
  structural_innovation_corr_prior_vector.clear();
  measurement_error_corr_prior_vector.clear();
  exo_names.clear();
  exo_det_names.clear();
  endo_names.clear();
  param_names.clear();
  params.clear();
  aux_vars.clear();
  predetermined_variables.clear();
  varobs.clear();
  lead_lag_incidence.clear();
  NNZDerivatives.clear();
}

string
DynareInfo::get_exo_name_by_index(int index) throw (ValueNotSetException)
{
  for (map<string, int >::iterator it = exo_names.begin();
       it != exo_names.end(); it++)
    if (it->second == index)
      return it->first;
  throw ValueNotSetException("get_exo_name_by_index" + index);
}

int
DynareInfo::get_exo_index_by_name(string name) throw (ValueNotSetException)
{
  map<string, int >::iterator it = exo_names.find(name);
  if (it != exo_names.end())
    return it->second;
  throw ValueNotSetException("get_exo_name_by_name" + name);
}

string
DynareInfo::get_exo_det_name_by_index(int index) throw (ValueNotSetException)
{
  for (map<string, int >::iterator it = exo_det_names.begin();
       it != exo_det_names.end(); it++)
    if (it->second == index)
      return it->first;
  throw ValueNotSetException("get_exo_det_name_by_index" + index);
}

int
DynareInfo::get_exo_det_index_by_name(string name) throw (ValueNotSetException)
{
  map<string, int >::iterator it = exo_det_names.find(name);
  if (it != exo_det_names.end())
    return it->second;
  throw ValueNotSetException("get_exo_det_name_by_name" + name);
}

string
DynareInfo::get_endo_name_by_index(int index) throw (ValueNotSetException)
{
  for (map<string, int >::iterator it = endo_names.begin();
       it != endo_names.end(); it++)
    if (it->second == index)
      return it->first;
  throw ValueNotSetException("get_endo_name_by_index" + index);
}

int
DynareInfo::get_endo_index_by_name(string name) throw (ValueNotSetException)
{
  map<string, int >::iterator it = endo_names.find(name);
  if (it != endo_names.end())
    return it->second;
  throw ValueNotSetException("get_endo_name_by_name" + name);
}

string
DynareInfo::get_param_name_by_index(int index) throw (ValueNotSetException)
{
  for (map<string, int >::iterator it = param_names.begin();
       it != param_names.end(); it++)
    if (it->second == index)
      return it->first;
  throw ValueNotSetException("get_param_name_by_index" + index);
}

int
DynareInfo::get_param_index_by_name(string name) throw (ValueNotSetException)
{
  map<string, int >::iterator it = param_names.find(name);
  if (it != param_names.end())
    return it->second;
  throw ValueNotSetException("get_param_name_by_name" + name);
}

double
DynareInfo::get_param_value_by_index(int index) throw (ValueNotSetException)
{
  return params[index];
  //  map<int, double >::iterator it = params.find(index);
  //  if (it != params.end())
  //    return it->second;
  //  throw ValueNotSetException("get_param_value_by_index" + index);
}

vector<int >
DynareInfo::get_lead_lag_incidence_for_endo_var_by_index(int index) throw (ValueNotSetException)
{
  if (index < lead_lag_incidence.size())
    return lead_lag_incidence.at(index);
  throw ValueNotSetException("get_lead_lag_incidence_for_endo_var_by_index" + index);
}

MarkovSwitching::MarkovSwitching(const int chain_arg,
                                 const int number_of_regimes_arg,
                                 const int number_of_lags_arg,
                                 const bool number_of_lags_was_passed_arg,
                                 const vector<int> parameters_arg,
                                 const vector<double> duration_arg,
                                 const restriction_map_t restriction_map_arg) :
  chain(chain_arg),
  number_of_regimes(number_of_regimes_arg),
  number_of_lags(number_of_lags_arg),
  number_of_lags_was_passed(number_of_lags_was_passed_arg),
  parameters(parameters_arg),
  duration(duration_arg),
  restriction_map(restriction_map_arg)
{
  assert(chain >= 1);
  assert(number_of_regimes > 0);
  if (number_of_lags_was_passed)
    assert(number_of_lags > 0);
  assert(!parameters.empty());
  assert(!duration.empty());
}

int
MarkovSwitching::get_number_of_lags() throw (ValueNotSetException)
{
  if (number_of_lags_has_val())
    return number_of_lags;
  throw ValueNotSetException("number_of_lags");
}

restriction_map_t
MarkovSwitching::get_restriction_map() throw (ValueNotSetException)
{
  if (restriction_map_has_val())
    return restriction_map;
  throw ValueNotSetException("restriction_map");
}

BasicModFilePrior::BasicModFilePrior(const int index_arg,
                                     const string shape_arg,
                                     const double mean_arg,
                                     const double mode_arg,
                                     const double stdev_arg,
                                     const double variance_arg,
                                     const vector <double> domain_arg) :
  index(index_arg),
  shape(shape_arg),
  mean(mean_arg),
  mode(mode_arg),
  stdev(stdev_arg),
  variance(variance_arg),
  domain(domain_arg)
{
  assert(index >= 0);
  assert(!shape.empty());
  if (stdev_has_val())
    assert(stdev >= 0);
  if (variance_has_val())
    assert(variance >= 0);
  if (domain_has_val())
    assert(domain.size() == 2);
}

double
BasicModFilePrior::get_mean() throw (ValueNotSetException)
{
  if (mean_has_val())
    return mean;
  throw ValueNotSetException("mean");
};

double
BasicModFilePrior::get_mode() throw (ValueNotSetException)
{
  if (mode_has_val())
    return mode;
  throw ValueNotSetException("mode");
};

double
BasicModFilePrior::get_stdev() throw (ValueNotSetException)
{
  if (stdev_has_val())
    return stdev;
  throw ValueNotSetException("stdev");
};

double
BasicModFilePrior::get_variance() throw (ValueNotSetException)
{
  if (variance_has_val())
    return variance;
  throw ValueNotSetException("variance");
};

vector<double>
BasicModFilePrior::get_domain() throw (ValueNotSetException)
{
  if (domain_has_val())
    return domain;
  throw ValueNotSetException("domain");
};

ModFilePrior::ModFilePrior(const int index_arg,
                           const string shape_arg,
                           const double mean_arg,
                           const double mode_arg,
                           const double stdev_arg,
                           const double variance_arg,
                           const vector <double> domain_arg) :
  BasicModFilePrior(index_arg,
                    shape_arg,
                    mean_arg,
                    mode_arg,
                    stdev_arg,
                    variance_arg,
                    domain_arg)
{
}

ModFileStructuralInnovationPrior::ModFileStructuralInnovationPrior(const int index_arg,
                                                                   const string shape_arg,
                                                                   const double mean_arg,
                                                                   const double stdev_arg,
                                                                   const double variance_arg,
                                                                   const double mode_arg,
                                                                   const vector <double> domain_arg) :
  BasicModFilePrior(index_arg,
                    shape_arg,
                    mean_arg,
                    mode_arg,
                    stdev_arg,
                    variance_arg,
                    domain_arg)
{
}

ModFileMeasurementErrorPrior::ModFileMeasurementErrorPrior(const int index_arg,
                                                           const string shape_arg,
                                                           const double mean_arg,
                                                           const double stdev_arg,
                                                           const double variance_arg,
                                                           const double mode_arg,
                                                           const vector <double> domain_arg) :
  BasicModFilePrior(index_arg,
                    shape_arg,
                    mean_arg,
                    mode_arg,
                    stdev_arg,
                    variance_arg,
                    domain_arg)
{
}

ModFileStructuralInnovationCorrPrior::ModFileStructuralInnovationCorrPrior(const int index1_arg,
                                                                           const int index2_arg,
                                                                           const string shape_arg,
                                                                           const double mean_arg,
                                                                           const double stdev_arg,
                                                                           const double variance_arg,
                                                                           const double mode_arg,
                                                                           const vector <double> domain_arg) :
  BasicModFilePrior(index1_arg,
                    shape_arg,
                    mean_arg,
                    mode_arg,
                    stdev_arg,
                    variance_arg,
                    domain_arg),
  index2(index2_arg)
{
  assert(index2 >= 0);
}

ModFileMeasurementErrorCorrPrior::ModFileMeasurementErrorCorrPrior(const int index1_arg,
                                                                   const int index2_arg,
                                                                   const string shape_arg,
                                                                   const double mean_arg,
                                                                   const double stdev_arg,
                                                                   const double variance_arg,
                                                                   const double mode_arg,
                                                                   const vector <double> domain_arg) :
  BasicModFilePrior(index1_arg,
                    shape_arg,
                    mean_arg,
                    mode_arg,
                    stdev_arg,
                    variance_arg,
                    domain_arg),
  index2(index2_arg)
{
  assert(index2 >= 0);
}

BasicModFileOption::BasicModFileOption(const int index_arg,
                                       const double init_arg) :
  index(index_arg),
  init(init_arg)
{
  assert(index >= 0);
  assert(!isnan(init));
}

ModFileOption::ModFileOption(const int index_arg, const double init_arg) :
  BasicModFileOption(index_arg, init_arg)
{
}

ModFileStructuralInnovationOption::ModFileStructuralInnovationOption(const int index_arg, const double init_arg) :
  BasicModFileOption(index_arg, init_arg)
{
}

ModFileMeasurementErrorOption::ModFileMeasurementErrorOption(const int index_arg, const double init_arg) :
  BasicModFileOption(index_arg, init_arg)
{
}

ModFileStructuralInnovationCorrOption::ModFileStructuralInnovationCorrOption(const int index1_arg, const int index2_arg, const double init_arg) :
  BasicModFileOption(index1_arg, init_arg),
  index2(index2_arg)
{
  assert(index2 >= 0);
}

ModFileMeasurementErrorCorrOption::ModFileMeasurementErrorCorrOption(const int index1_arg, const int index2_arg, const double init_arg) :
  BasicModFileOption(index1_arg, init_arg),
  index2(index2_arg)
{
  assert(index2 >= 0);
}
