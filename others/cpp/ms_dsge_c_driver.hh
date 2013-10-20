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

#ifndef _MS_DSGE_C_DRIVER_HH
#define _MS_DSGE_C_DRIVER_HH

#include <cstdlib>
#include <vector>
#include <string>
#include <map>
#include <limits>

using namespace std;

struct aux_vars_t {
  int endo_index, type, orig_index, orig_lead_lag;
} ;

typedef map<pair<int, int >, double> restriction_map_t ;

class ValueNotSetException
{
public:
  string name;
  ValueNotSetException(const string &name_arg) : name(name_arg)
  {
  }
};

class MarkovSwitching
{
private:
  const int chain, number_of_regimes, number_of_lags;
  const bool number_of_lags_was_passed;
  const vector<int> parameters;
  const vector<double> duration;
  const restriction_map_t restriction_map;
public:
  MarkovSwitching(const int chain_arg,
                  const int number_of_regimes_arg,
                  const int number_of_lags_arg,
                  const bool number_of_lags_was_passed_arg,
                  const vector<int> parameters_arg,
                  const vector<double> duration_arg,
                  const restriction_map_t restriction_map_arg);

  inline bool number_of_lags_has_val() const { return number_of_lags_was_passed; };
  inline bool restriction_map_has_val() const { return !restriction_map.empty(); };

  inline int get_chain() const { return chain; };
  inline int get_number_of_regimes() const { return number_of_regimes; };
  int get_number_of_lags() throw (ValueNotSetException);
  inline vector<int> get_parameters() { return parameters; };
  inline vector<double> get_duration() { return duration; };
  restriction_map_t get_restriction_map() throw (ValueNotSetException);
};

class BasicModFilePrior
{
private:
  inline bool isnan(double x) const { return (x!=x); };
protected:
  const int index;
  const string shape;
  const double mean, mode, stdev, variance;
  const vector <double> domain;

  BasicModFilePrior(const int index_arg,
                    const string shape_arg,
                    const double mean_arg,
                    const double mode_arg,
                    const double stdev_arg,
                    const double variance_arg,
                    const vector <double> domain_arg);
public:
  inline bool mean_has_val() const { return !isnan(mean); };
  inline bool mode_has_val() const { return !isnan(mode); };
  inline bool stdev_has_val() const { return !isnan(stdev); };
  inline bool variance_has_val() const { return !isnan(variance); };
  inline bool domain_has_val() const { return (domain.size() == 2); };

  inline int get_index() const { return index; };
  inline string get_shape() const { return shape; };
  double get_mean() throw (ValueNotSetException);
  double get_mode() throw (ValueNotSetException);
  double get_stdev() throw (ValueNotSetException);
  double get_variance() throw (ValueNotSetException);
  vector<double> get_domain() throw (ValueNotSetException);
};

class ModFilePrior : public BasicModFilePrior
{
public:
  ModFilePrior(const int index_arg,
               const string shape_arg,
               const double mean_arg,
               const double mode_arg,
               const double stdev_arg,
               const double variance_arg,
               const vector <double> domain_arg);
};

class ModFileStructuralInnovationPrior: public BasicModFilePrior
{
public:
  ModFileStructuralInnovationPrior(const int index_arg,
                                   const string shape_arg,
                                   const double mean_arg,
                                   const double mode_arg,
                                   const double stdev_arg,
                                   const double variance_arg,
                                   const vector <double> domain_arg);
};

class ModFileMeasurementErrorPrior : public BasicModFilePrior
{
public:
  ModFileMeasurementErrorPrior(const int index_arg,
                               const string shape_arg,
                               const double mean_arg,
                               const double stdev_arg,
                               const double variance_arg,
                               const double mode_arg,
                               const vector <double> domain_arg);
};

class ModFileStructuralInnovationCorrPrior : public BasicModFilePrior
{
private:
  const int index2;
public:
  ModFileStructuralInnovationCorrPrior(const int index1_arg,
                                       const int index2_arg,
                                       const string shape_arg,
                                       const double mean_arg,
                                       const double stdev_arg,
                                       const double variance_arg,
                                       const double mode_arg,
                                       const vector <double> domain_arg);
};

class ModFileMeasurementErrorCorrPrior : public BasicModFilePrior
{
private:
  const int index2;
public:
  ModFileMeasurementErrorCorrPrior(const int index1_arg,
                                   const int index2_arg,
                                   const string shape_arg,
                                   const double mean_arg,
                                   const double stdev_arg,
                                   const double variance_arg,
                                   const double mode_arg,
                                   const vector <double> domain_arg);
};

class BasicModFileOption
{
private:
  inline bool isnan(double x) const { return (x!=x); };
protected:
  const int index;
  const double init;

  BasicModFileOption(const int index_arg,
                     const double init_arg);
public:
  inline int get_index() const { return index; };
  inline double get_init() const { return init; };
};

class ModFileOption : public BasicModFileOption
{
public:
  ModFileOption(const int index_arg,
                const double init_arg);
};

class ModFileStructuralInnovationOption: public BasicModFileOption
{
public:
  ModFileStructuralInnovationOption(const int index_arg,
                                    const double init_arg);
};

class ModFileMeasurementErrorOption : public BasicModFileOption
{
public:
  ModFileMeasurementErrorOption(const int index_arg,
                                const double init_arg);
};

class ModFileStructuralInnovationCorrOption : public BasicModFileOption
{
private:
  const int index2;
public:
  ModFileStructuralInnovationCorrOption(const int index1_arg,
                                        const int index2_arg,
                                        const double init_arg);
};

class ModFileMeasurementErrorCorrOption : public BasicModFileOption
{
private:
  const int index2;
public:
  ModFileMeasurementErrorCorrOption(const int index1_arg,
                                    const int index2_arg,
                                    const double init_arg);
};

class MsDsgeInfo
{
private:
  vector<MarkovSwitching *> markov_switching_vector;
  vector<ModFilePrior *> prior_vector;
  vector<ModFileStructuralInnovationPrior *> structural_innovation_prior_vector;
  vector<ModFileMeasurementErrorPrior *> measurement_error_prior_vector;
  vector<ModFileStructuralInnovationCorrPrior *> structural_innovation_corr_prior_vector;
  vector<ModFileMeasurementErrorCorrPrior *> measurement_error_corr_prior_vector;
  vector<ModFileOption *> option_vector;
  vector<ModFileStructuralInnovationOption *> structural_innovation_option_vector;
  vector<ModFileMeasurementErrorOption *> measurement_error_option_vector;
  vector<ModFileStructuralInnovationCorrOption *> structural_innovation_corr_option_vector;
  vector<ModFileMeasurementErrorCorrOption *> measurement_error_corr_option_vector;
  map<string, int > exo_names, exo_det_names, endo_names, param_names;
  map<int, double > params;
  vector<aux_vars_t> aux_vars;
  vector<int> predetermined_variables;
  vector<int> varobs;
  vector<vector<int > >lead_lag_incidence;
  vector<double> NNZDerivatives;
public:
  MsDsgeInfo(map<string, int > exo_names_arg,
             map<string, int > exo_det_names_arg,
             map<string, int > endo_names_arg,
             map<string, int > param_names_arg,
             map<int, double > params_arg,
             vector<aux_vars_t> aux_vars_arg,
             vector<int> predetermined_variables_arg,
             vector<int> varobs_arg,
             vector<vector<int > > lead_lag_incidence_arg,
             vector<double> NNZDerivatives_arg);
  ~MsDsgeInfo();

  inline void addMarkovSwitching(MarkovSwitching *ms) { markov_switching_vector.push_back(ms); };

  inline void addPrior(ModFilePrior *p) { prior_vector.push_back(p); };
  inline void addStructuralInnovationPrior(ModFileStructuralInnovationPrior *sip) { structural_innovation_prior_vector.push_back(sip); };
  inline void addMeasurementErrorPrior(ModFileMeasurementErrorPrior *mep) { measurement_error_prior_vector.push_back(mep); };
  inline void addStructuralInnovationCorrPrior(ModFileStructuralInnovationCorrPrior *sicp) { structural_innovation_corr_prior_vector.push_back(sicp); };
  inline void addMeasurementErrorCorrPrior(ModFileMeasurementErrorCorrPrior *mecp) { measurement_error_corr_prior_vector.push_back(mecp); };

  inline void addOption(ModFileOption *o) { option_vector.push_back(o); };
  inline void addStructuralInnovationOption(ModFileStructuralInnovationOption *sio) { structural_innovation_option_vector.push_back(sio); };
  inline void addMeasurementErrorOption(ModFileMeasurementErrorOption *meo) { measurement_error_option_vector.push_back(meo); };
  inline void addStructuralInnovationCorrOption(ModFileStructuralInnovationCorrOption *sico) { structural_innovation_corr_option_vector.push_back(sico); };
  inline void addMeasurementErrorCorrOption(ModFileMeasurementErrorCorrOption *meco) { measurement_error_corr_option_vector.push_back(meco); };

  inline bool markov_switching_has_val() { return !markov_switching_vector.empty(); };
  inline bool prior_has_val() { return !prior_vector.empty(); };
  inline bool structural_innovation_prior_has_val() { return !structural_innovation_prior_vector.empty(); };
  inline bool measurement_error_prior_has_val() { return !measurement_error_prior_vector.empty(); };
  inline bool structural_innovation_corr_prior_has_val() { return !structural_innovation_corr_prior_vector.empty(); };
  inline bool measurement_error_corr_prior_has_val() { return !measurement_error_corr_prior_vector.empty(); };

  inline bool option_has_val() { return !option_vector.empty(); };
  inline bool structural_innovation_option_has_val() { return !structural_innovation_option_vector.empty(); };
  inline bool measurement_error_option_has_val() { return !measurement_error_option_vector.empty(); };
  inline bool structural_innovation_corr_option_has_val() { return !structural_innovation_corr_option_vector.empty(); };
  inline bool measurement_error_corr_option_has_val() { return !measurement_error_corr_option_vector.empty(); };

  inline vector<MarkovSwitching *>get_markov_switching() { return markov_switching_vector; };
  inline vector<ModFilePrior *> get_prior() { return prior_vector; };
  inline vector<ModFileStructuralInnovationPrior *> get_structural_innovation_prior() { return structural_innovation_prior_vector; };
  inline vector<ModFileMeasurementErrorPrior *> get_measurement_error_prior() { return measurement_error_prior_vector; };
  inline vector<ModFileStructuralInnovationCorrPrior *> get_structural_innovation_corr_prior() { return structural_innovation_corr_prior_vector; };
  inline vector<ModFileMeasurementErrorCorrPrior *> get_measurement_error_corr_prior() { return measurement_error_corr_prior_vector; };

  inline vector<ModFileOption *> get_option() { return option_vector; };
  inline vector<ModFileStructuralInnovationOption *> get_structural_innovation_option() { return structural_innovation_option_vector; };
  inline vector<ModFileMeasurementErrorOption *> get_measurement_error_option() { return measurement_error_option_vector; };
  inline vector<ModFileStructuralInnovationCorrOption *> get_structural_innovation_corr_option() { return structural_innovation_corr_option_vector; };
  inline vector<ModFileMeasurementErrorCorrOption *> get_measurement_error_corr_option() { return measurement_error_corr_option_vector; };

  inline map<string, int > get_exo_names() { return exo_names; };
  inline map<string, int > get_exo_det_names() { return exo_det_names; };
  inline map<string, int > get_endo_names() { return endo_names; };
  inline map<string, int > get_param_names() { return param_names; };
  inline map<int, double > get_params() { return params; };
  inline vector <aux_vars_t> get_aux_vars() { return aux_vars; };
  inline vector <int> get_predetermined_variables() { return predetermined_variables; };
  inline vector <int> get_varobs() { return varobs; };
  inline vector<vector<int > > get_lead_lag_incidence() { return lead_lag_incidence; };
  inline vector<double> get_NNZDerivatives() { return NNZDerivatives; };

  string get_exo_name_by_index(int index) throw (ValueNotSetException);
  int get_exo_index_by_name(string name) throw (ValueNotSetException);
  string get_exo_det_name_by_index(int index) throw (ValueNotSetException);
  int get_exo_det_index_by_name(string name) throw (ValueNotSetException);
  string get_endo_name_by_index(int index) throw (ValueNotSetException);
  int get_endo_index_by_name(string name) throw (ValueNotSetException);
  string get_param_name_by_index(int index) throw (ValueNotSetException);
  int get_param_index_by_name(string name) throw (ValueNotSetException);
  double get_param_value_by_index(int index) throw (ValueNotSetException);
  vector<int >get_lead_lag_incidence_for_endo_var_by_index(int index) throw (ValueNotSetException);
};

#endif // ! _MS_DSGE_C_DRIVER_HH
