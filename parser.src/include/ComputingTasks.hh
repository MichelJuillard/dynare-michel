#ifndef _COMPUTINGTASKS_HH
#define _COMPUTINGTASKS_HH

#include <ostream>

#include "TmpSymbolTable.hh"
#include "SymbolTable.hh"
#include "Statement.hh"

class SteadyStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  SteadyStatement(const OptionsList &options_list_arg);
  virtual void writeOutput(ostream &output) const;
};

class CheckStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  CheckStatement(const OptionsList &options_list_arg);
  virtual void writeOutput(ostream &output) const;
};

class SimulStatement : public Statement
{
private:
  const OptionsList options_list;
public:
  SimulStatement(const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output) const;
};

class StochSimulStatement : public Statement
{
private:
  const TmpSymbolTable tmp_symbol_table;
  const OptionsList options_list;
public:
  StochSimulStatement(const TmpSymbolTable &tmp_symbol_table_arg,
                      const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output) const;
};

class RplotStatement : public Statement
{
private:
  const TmpSymbolTable tmp_symbol_table;
  const OptionsList options_list;
public:
  RplotStatement(const TmpSymbolTable &tmp_symbol_table_arg,
                 const OptionsList &options_list_arg);
  virtual void writeOutput(ostream &output) const;
};

class UnitRootVarsStatement : public Statement
{
private:
  const TmpSymbolTable tmp_symbol_table;
public:
  UnitRootVarsStatement(const TmpSymbolTable &tmp_symbol_table_arg);
  virtual void writeOutput(ostream &output) const;
};

class PeriodsStatement : public Statement
{
private:
  const int periods;
public:
  PeriodsStatement(int periods_arg);
  virtual void writeOutput(ostream &output) const;
};

class DsampleStatement : public Statement
{
private:
  const int val1, val2;
public:
  DsampleStatement(int val1_arg);
  DsampleStatement(int val1_arg, int val2_arg);
  virtual void writeOutput(ostream &output) const;
};

class EstimationStatement : public Statement
{
private:
  const TmpSymbolTable tmp_symbol_table;
  const OptionsList options_list;
public:
  EstimationStatement(const TmpSymbolTable &tmp_symbol_table_arg,
                      const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output) const;
};

class PriorAnalysisStatement : public Statement
{
private:
  const TmpSymbolTable tmp_symbol_table;
  const OptionsList options_list;
public:
  PriorAnalysisStatement(const TmpSymbolTable &tmp_symbol_table_arg,
                         const OptionsList &options_list_arg);
  virtual void writeOutput(ostream &output) const;
};

class PosteriorAnalysisStatement : public Statement
{
private:
  const TmpSymbolTable tmp_symbol_table;
  const OptionsList options_list;
public:
  PosteriorAnalysisStatement(const TmpSymbolTable &tmp_symbol_table_arg,
                             const OptionsList &options_list_arg);
  virtual void writeOutput(ostream &output) const;
};

class VarobsStatement : public Statement
{
private:
  const TmpSymbolTable tmp_symbol_table;
public:
  VarobsStatement(const TmpSymbolTable &tmp_symbol_table_arg);
  virtual void writeOutput(ostream &output) const;
};

class ObservationTrendsStatement : public Statement
{
public:
  typedef map<string, string, less<string> > trend_elements_type;
private:
  const trend_elements_type trend_elements;
  const SymbolTable &symbol_table;
public:
  ObservationTrendsStatement(const trend_elements_type &trend_elements_arg,
                             const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output) const;
};

class OsrParamsStatement : public Statement
{
private:
  const TmpSymbolTable tmp_symbol_table;
public:
  OsrParamsStatement(const TmpSymbolTable &tmp_symbol_table_arg);
  virtual void writeOutput(ostream &output) const;
};

class OsrStatement : public Statement
{
private:
  const TmpSymbolTable tmp_symbol_table;
  const OptionsList options_list;
public:
  OsrStatement(const TmpSymbolTable &tmp_symbol_table_arg,
               const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output) const;
};

class OlrStatement : public Statement
{
private:
  const TmpSymbolTable tmp_symbol_table;
  const OptionsList options_list;
public:
  OlrStatement(const TmpSymbolTable &tmp_symbol_table_arg,
               const OptionsList &options_list_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct);
  virtual void writeOutput(ostream &output) const;
};

class OlrInstStatement : public Statement
{
private:
  const TmpSymbolTable tmp_symbol_table;
public:
  OlrInstStatement(const TmpSymbolTable &tmp_symbol_table_arg);
  virtual void writeOutput(ostream &output) const;
};

class DynaTypeStatement : public Statement
{
private:
  const TmpSymbolTable tmp_symbol_table;
  const string filename;
  const string ext;
public:
  DynaTypeStatement(const TmpSymbolTable &tmp_symbol_table_arg,
                    const string &filename_arg, const string &ext_arg);
  virtual void writeOutput(ostream &output) const;
};

class DynaSaveStatement : public Statement
{
private:
  const TmpSymbolTable tmp_symbol_table;
  const string filename;
  const string ext;
public:
  DynaSaveStatement(const TmpSymbolTable &tmp_symbol_table_arg,
                    const string &filename_arg, const string &ext_arg);
  virtual void writeOutput(ostream &output) const;
};

class ModelComparisonStatement : public Statement
{
public:
  typedef map<string, string, less<string> > filename_list_type;
private:
  filename_list_type filename_list;
  OptionsList options_list;
public:
  ModelComparisonStatement(const filename_list_type &filename_list_arg,
                           const OptionsList &options_list_arg);
  virtual void writeOutput(ostream &output) const;
};

/*!
  \class EstimationParams
  \brief EstimationParams
*/
struct EstimationParams
{
  int type;
  std::string name;
  std::string name2;
  std::string init_val;
  std::string prior;
  std::string low_bound;
  std::string up_bound;
  std::string mean;
  std::string std;
  std::string p3;
  std::string p4;
  std::string jscale;

  EstimationParams()
  {
    clear();
  }
  void clear()
  {
    type = 0;
    name = "";
    name2 = "";
    init_val = "NaN";
    prior = "NaN";
    low_bound = "-Inf";
    up_bound = "Inf";
    mean = "NaN";
    std = "NaN";
    p3 = "NaN";
    p4 = "NaN";
    jscale = "NaN";
  }
};

class EstimatedParamsStatement : public Statement
{
private:
  const vector<EstimationParams> estim_params_list;
  const SymbolTable &symbol_table;
public:
  EstimatedParamsStatement(const vector<EstimationParams> &estim_params_list_arg,
                           const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output) const;
};

class EstimatedParamsInitStatement : public Statement
{
private:
  const vector<EstimationParams> estim_params_list;
  const SymbolTable &symbol_table;
public:
  EstimatedParamsInitStatement(const vector<EstimationParams> &estim_params_list_arg,
                               const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output) const;
};

class EstimatedParamsBoundsStatement : public Statement
{
private:
  const vector<EstimationParams> estim_params_list;
  const SymbolTable &symbol_table;
public:
  EstimatedParamsBoundsStatement(const vector<EstimationParams> &estim_params_list_arg,
                                 const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output) const;
};

class OptimWeightsStatement : public Statement
{
public:
  typedef map<string, string, less<string> > var_weights_type;
  typedef map<pair<string, string>, string, less<pair<string, string> > > covar_weights_type;
private:
  const var_weights_type var_weights;
  const covar_weights_type covar_weights;
  const SymbolTable &symbol_table;
public:
  OptimWeightsStatement(const var_weights_type &var_weights_arg,
                        const covar_weights_type &covar_weights_arg,
                        const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output) const;
};

class CalibStatement : public Statement
{
private:
  const int covar;
public:
  CalibStatement(int covar_arg);
  virtual void writeOutput(ostream &output) const;
};

class CalibVarStatement : public Statement
{
public:
  //! Maps a variable to a pair (weight, expression)
  typedef map<string, pair<string, string>, less<string> > calib_var_type;
  //! Maps a pair of variables to a pair (weight, expression)
  typedef map<pair<string, string>, pair<string, string>, less<pair<string, string> > > calib_covar_type;
  //! Maps a pair (variable, autocorr) to a pair (weight, expression)
  typedef map<pair<string, int>, pair<string, string>, less<pair<string, int> > > calib_ac_type;
private:
  const calib_var_type calib_var;
  const calib_covar_type calib_covar;
  const calib_ac_type calib_ac;
  const SymbolTable &symbol_table;
public:
  CalibVarStatement(const calib_var_type &calib_var_arg,
                    const calib_covar_type &calib_covar_arg,
                    const calib_ac_type &calib_ac_arg,
                    const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output) const;
};

#endif
