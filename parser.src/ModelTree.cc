#include <iostream>
#include <fstream>
#include <sstream>

#include "ModelTree.hh"
#include "Interface.hh"

ModelTree::ModelTree(SymbolTable &symbol_table_arg,
                     NumericalConstants &num_constants_arg) :
  DataTree(symbol_table_arg, num_constants_arg),
  computeJacobian(false),
  computeJacobianExo(false),
  computeHessian(false),
  computeStaticHessian(false),
  computeThirdDerivatives(false)
{
}

void
ModelTree::derive(int order)
{
  cout << "Processing derivation ..." << endl;

  cout << "  Processing Order 1... ";
  for(int var = 0; var < variable_table.size(); var++)
    for(int eq = 0; eq < (int) equations.size(); eq++)
      {
        NodeID d1 = equations[eq]->getDerivative(var);
        if (d1 == Zero)
          continue;
        first_derivatives[make_pair(eq, var)] = d1;
      }
  cout << "done" << endl;

  if (order >= 2)
    {
      cout << "  Processing Order 2... ";
      for(first_derivatives_type::const_iterator it = first_derivatives.begin();
          it != first_derivatives.end(); it++)
        {
          int eq = it->first.first;
          int var1 = it->first.second;
          NodeID d1 = it->second;
      
          // Store only second derivatives with var2 <= var1
          for(int var2 = 0; var2 <= var1; var2++)
            {
              NodeID d2 = d1->getDerivative(var2);
              if (d2 == Zero)
                continue;
              second_derivatives[make_pair(eq, make_pair(var1, var2))] = d2;
            }
        }
      cout << "done" << endl;
    }

  if (order >= 3)
    {
      cout << "  Processing Order 3... ";
      for(second_derivatives_type::const_iterator it = second_derivatives.begin();
          it != second_derivatives.end(); it++)
        {
          int eq = it->first.first;

          int var1 = it->first.second.first;
          int var2 = it->first.second.second;
          // By construction, var2 <= var1

          NodeID d2 = it->second;

          // Store only third derivatives such that var3 <= var2 <= var1
          for(int var3 = 0; var3 <= var2; var3++)
            {
              NodeID d3 = d2->getDerivative(var3);
              if (d3 == Zero)
                continue;
              third_derivatives[make_pair(eq, make_pair(var1, make_pair(var2, var3)))] = d3;
            }
        }
      cout << "done" << endl;
    }
}

void
ModelTree::computeTemporaryTerms(int order)
{
  map<NodeID, int> reference_count;
  temporary_terms.clear();

  for(vector<BinaryOpNode *>::iterator it = equations.begin();
      it != equations.end(); it++)
    (*it)->computeTemporaryTerms(reference_count, temporary_terms);

  for(first_derivatives_type::iterator it = first_derivatives.begin();
      it != first_derivatives.end(); it++)
    it->second->computeTemporaryTerms(reference_count, temporary_terms);

  if (order >= 2)
    for(second_derivatives_type::iterator it = second_derivatives.begin();
        it != second_derivatives.end(); it++)
      it->second->computeTemporaryTerms(reference_count, temporary_terms);

  if (order >= 3)
    for(third_derivatives_type::iterator it = third_derivatives.begin();
        it != third_derivatives.end(); it++)
      it->second->computeTemporaryTerms(reference_count, temporary_terms);
}

void
ModelTree::writeTemporaryTerms(ostream &output, bool is_dynamic) const
{
  // A copy of temporary terms
  temporary_terms_type tt2;

  for(temporary_terms_type::const_iterator it = temporary_terms.begin();
      it != temporary_terms.end(); it++)
    {
      (*it)->writeOutput(output, is_dynamic, temporary_terms);
      output << " = ";

      (*it)->writeOutput(output, is_dynamic, tt2);

      // Insert current node into tt2
      tt2.insert(*it);

      output << ";" << endl;
    }
}

void
ModelTree::writeLocalParameters(ostream &output, bool is_dynamic) const
{
  for(map<int, NodeID>::const_iterator it = local_parameters_table.begin();
      it != local_parameters_table.end(); it++)
    {
      int id = it->first;
      NodeID value = it->second;
      output << symbol_table.getNameByID(eLocalParameter, id) << " = ";
      value->writeOutput(output, is_dynamic, temporary_terms);
      output << ";" << endl;
    }
}

void
ModelTree::writeModelEquations(ostream &output, bool is_dynamic) const
{
  for(int eq = 0; eq < (int) equations.size(); eq++)
    {
      BinaryOpNode *eq_node = equations[eq];

      NodeID lhs = eq_node->arg1;
      output << "lhs =";
      lhs->writeOutput(output, is_dynamic, temporary_terms);
      output << ";" << endl;

      NodeID rhs = eq_node->arg2;
      output << "rhs =";
      rhs->writeOutput(output, is_dynamic, temporary_terms);
      output << ";" << endl;

      output << "residual" << lpar << eq + 1 << rpar << "= lhs-rhs;" << endl;
    }
}

void
ModelTree::writeStaticMFile(const string &static_basename) const
{
  string filename = static_basename + interfaces::function_file_extension();

  ofstream mStaticModelFile;
  mStaticModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mStaticModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(-1);
    }
  // Writing comments and function definition command
  mStaticModelFile << "function [residual, g1, g2] = " << static_basename << "( y, x )\n";
  mStaticModelFile << interfaces::comment()+"\n"+interfaces::comment();
  mStaticModelFile << "Status : Computes static model for Dynare\n" << interfaces::comment() << "\n";
  mStaticModelFile << interfaces::comment();
  mStaticModelFile << "Warning : this file is generated automatically by Dynare\n";
  mStaticModelFile << interfaces::comment();
  mStaticModelFile << "  from model file (.mod)\n\n";

  writeStaticModel(mStaticModelFile);

  interfaces::function_close();
  mStaticModelFile.close();
}


void
ModelTree::writeDynamicMFile(const string &dynamic_basename) const
{
  string filename = dynamic_basename + interfaces::function_file_extension();

  ofstream mDynamicModelFile;
  mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(-1);
    }
  mDynamicModelFile << "function [residual, g1, g2, g3] = " << dynamic_basename << "(y, x)\n";
  mDynamicModelFile << interfaces::comment()+"\n"+interfaces::comment();
  mDynamicModelFile << "Status : Computes dynamic model for Dynare\n" << interfaces::comment() << "\n";
  mDynamicModelFile << interfaces::comment();
  mDynamicModelFile << "Warning : this file is generated automatically by Dynare\n";
  mDynamicModelFile << interfaces::comment();
  mDynamicModelFile << "  from model file (.mod)\n\n";

  writeDynamicModel(mDynamicModelFile);  

  interfaces::function_close();
  mDynamicModelFile.close();
}

void
ModelTree::writeStaticCFile(const string &static_basename) const
{
  string filename = static_basename + ".c";

  ofstream mStaticModelFile;
  mStaticModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mStaticModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(-1);
    }
  mStaticModelFile << "/*\n";
  mStaticModelFile << " * " << filename << " : Computes static model for Dynare\n";
  mStaticModelFile << " * Warning : this file is generated automatically by Dynare\n";
  mStaticModelFile << " *           from model file (.mod)\n\n";
  mStaticModelFile << " */\n";
  mStaticModelFile << "#include <math.h>\n";
  mStaticModelFile << "#include \"mex.h\"\n";
  // A flobal variable for model parameters
  mStaticModelFile << "double *params;\n";

  // Writing the function Static
  writeStaticModel(mStaticModelFile);

  // Writing the gateway routine
  mStaticModelFile << "/* The gateway routine */\n";
  mStaticModelFile << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])\n";
  mStaticModelFile << "{\n";
  mStaticModelFile << "  double *y, *x;\n";
  mStaticModelFile << "  double *residual, *g1;\n";
  mStaticModelFile << "  mxArray *M_;\n";
  mStaticModelFile << "\n";
  mStaticModelFile << "  /* Create a pointer to the input matrix y. */\n";
  mStaticModelFile << "  y = mxGetPr(prhs[0]);\n";
  mStaticModelFile << "\n";
  mStaticModelFile << "  /* Create a pointer to the input matrix x. */\n";
  mStaticModelFile << "  x = mxGetPr(prhs[1]);\n";
  mStaticModelFile << "\n";

  mStaticModelFile << "  residual = NULL;\n";
  mStaticModelFile << "  if (nlhs >= 1)\n";
  mStaticModelFile << "  {\n";
  mStaticModelFile << "      /* Set the output pointer to the output matrix residual. */\n";
  mStaticModelFile << "      plhs[0] = mxCreateDoubleMatrix(" << equations.size() << ",1, mxREAL);\n";
  mStaticModelFile << "     /* Create a C pointer to a copy of the output matrix residual. */\n";
  mStaticModelFile << "     residual = mxGetPr(plhs[0]);\n";
  mStaticModelFile << "  }\n\n";
  mStaticModelFile << "  g1 = NULL;\n";
  mStaticModelFile << "  if (nlhs >= 2)\n";
  mStaticModelFile << "  {\n";
  mStaticModelFile << "      /* Set the output pointer to the output matrix g1. */\n";
  mStaticModelFile << "      plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << symbol_table.endo_nbr << ", mxREAL);\n";
  mStaticModelFile << "      /* Create a C pointer to a copy of the output matrix g1. */\n";
  mStaticModelFile << "      g1 = mxGetPr(plhs[1]);\n";
  mStaticModelFile << "  }\n\n";
  mStaticModelFile << "  /* Gets model parameters from global workspace of Matlab */\n";
  mStaticModelFile << "  M_ = mexGetVariable(\"global\",\"M_\");\n";
  mStaticModelFile << "  if (M_ == NULL ){\n";
  mStaticModelFile << "	    mexPrintf(\"Global variable not found : \");\n";
  mStaticModelFile << "	    mexErrMsgTxt(\"M_ \\n\");\n";
  mStaticModelFile << "  }\n";
  mStaticModelFile << "  params = mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,\"params\")));\n";
  mStaticModelFile << "  /* Call the C Static. */\n";
  mStaticModelFile << "  Static(y, x, residual, g1);\n";
  mStaticModelFile << "}\n";
  mStaticModelFile.close();
}

 
void
ModelTree::writeDynamicCFile(const string &dynamic_basename) const
{
  string filename = dynamic_basename + ".c";

  ofstream mDynamicModelFile;
  mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(-1);
    }
  mDynamicModelFile << "/*\n";
  mDynamicModelFile << " * " << filename << " : Computes dynamic model for Dynare\n";
  mDynamicModelFile  << " *\n";
  mDynamicModelFile << " * Warning : this file is generated automatically by Dynare\n";
  mDynamicModelFile << " *           from model file (.mod)\n\n";
  mDynamicModelFile << " */\n";
  mDynamicModelFile << "#include <math.h>\n";
  mDynamicModelFile << "#include \"mex.h\"\n";
  // A flobal variable for model parameters
  mDynamicModelFile << "double *params;\n";
  // A global variable for it_
  mDynamicModelFile << "int it_;\n";
  mDynamicModelFile << "int nb_row_x;\n";

  // Writing the function body
  writeDynamicModel(mDynamicModelFile);  

  // Writing the gateway routine
  mDynamicModelFile << "/* The gateway routine */\n";
  mDynamicModelFile << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])\n";
  mDynamicModelFile << "{\n";
  mDynamicModelFile << "  double *y, *x;\n";
  mDynamicModelFile << "  double *residual, *g1, *g2;\n";
  mDynamicModelFile << "  mxArray *M_;\n";
  mDynamicModelFile << "\n";
  mDynamicModelFile << "  /* Create a pointer to the input matrix y. */\n";
  mDynamicModelFile << "  y = mxGetPr(prhs[0]);\n";
  mDynamicModelFile << "\n";
  mDynamicModelFile << "  /* Create a pointer to the input matrix x. */\n";
  mDynamicModelFile << "  x = mxGetPr(prhs[1]);\n";
  mDynamicModelFile << "  /* Gets number of rows of matrix x. */\n";
  mDynamicModelFile << "  nb_row_x = mxGetM(prhs[1]);\n";
  mDynamicModelFile << "\n";
  mDynamicModelFile << "  residual = NULL;\n";
  mDynamicModelFile << "  if (nlhs >= 1)\n";
  mDynamicModelFile << "  {\n";
  mDynamicModelFile << "     /* Set the output pointer to the output matrix residual. */\n";
  mDynamicModelFile << "     plhs[0] = mxCreateDoubleMatrix(" << equations.size() << ",1, mxREAL);\n";
  mDynamicModelFile << "     /* Create a C pointer to a copy of the output matrix residual. */\n";
  mDynamicModelFile << "     residual = mxGetPr(plhs[0]);\n";
  mDynamicModelFile << "  }\n\n";
  mDynamicModelFile << "  g1 = NULL;\n";
  mDynamicModelFile << "  if (nlhs >= 2)\n";
  mDynamicModelFile << "  {\n";
  mDynamicModelFile << "     /* Set the output pointer to the output matrix g1. */\n";
  if (computeJacobianExo)
    mDynamicModelFile << "     plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << variable_table.get_dyn_var_nbr() << ", mxREAL);\n";
  else if (computeJacobian)
    mDynamicModelFile << "     plhs[1] = mxCreateDoubleMatrix(" << equations.size() << ", " << variable_table.var_endo_nbr << ", mxREAL);\n";
  mDynamicModelFile << "     /* Create a C pointer to a copy of the output matrix g1. */\n";
  mDynamicModelFile << "     g1 = mxGetPr(plhs[1]);\n";
  mDynamicModelFile << "  }\n\n";
  mDynamicModelFile << "  g2 = NULL;\n";
  mDynamicModelFile << " if (nlhs >= 3)\n";
  mDynamicModelFile << "  {\n";
  mDynamicModelFile << "     /* Set the output pointer to the output matrix g2. */\n";
  mDynamicModelFile << "     plhs[2] = mxCreateDoubleMatrix(" << equations.size() << ", " << variable_table.get_dyn_var_nbr()*variable_table.get_dyn_var_nbr() << ", mxREAL);\n";
  mDynamicModelFile << "     /* Create a C pointer to a copy of the output matrix g1. */\n";
  mDynamicModelFile << "     g2 = mxGetPr(plhs[2]);\n";
  mDynamicModelFile << "  }\n\n";
  mDynamicModelFile << "  /* Gets model parameters from global workspace of Matlab */\n";
  mDynamicModelFile << "  M_ = mexGetVariable(\"global\",\"M_\");\n";
  mDynamicModelFile << "  if (M_ == NULL )\n";
  mDynamicModelFile << "  {\n";
  mDynamicModelFile << "	    mexPrintf(\"Global variable not found : \");\n";
  mDynamicModelFile << "	    mexErrMsgTxt(\"M_ \\n\");\n";
  mDynamicModelFile << "  }\n";
  mDynamicModelFile << "  params = mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,\"params\")));\n";
  mDynamicModelFile << "  /* Gets it_ from global workspace of Matlab */\n";
  mDynamicModelFile << "  it_ = (int) floor(mxGetScalar(mexGetVariable(\"global\", \"it_\")))-1;\n";
  mDynamicModelFile << "  /* Call the C subroutines. */\n";
  mDynamicModelFile << "  Dynamic(y, x, residual, g1, g2);\n";
  mDynamicModelFile << "}\n";
  mDynamicModelFile.close();
}

void
ModelTree::writeStaticModel(ostream &StaticOutput) const
{
  ostringstream model_output;    // Used for storing model equations
  ostringstream jacobian_output; // Used for storing jacobian equations
  ostringstream hessian_output;
  ostringstream lsymetric;       // For symmetric elements in hessian

  writeTemporaryTerms(model_output, false);

  writeLocalParameters(model_output, false);

  writeModelEquations(model_output, false);

  // Write Jacobian w.r. to endogenous only
  for(first_derivatives_type::const_iterator it = first_derivatives.begin();
      it != first_derivatives.end(); it++)
    {
      int eq = it->first.first;
      int var = it->first.second;
      NodeID d1 = it->second;

      if (variable_table.getType(var) == eEndogenous)
        {
          ostringstream g1;
          g1 << "  g1" << lpar << eq + 1 << ", " << variable_table.getSymbolID(var) + 1 << rpar;
          
          jacobian_output << g1.str() << "=" << g1.str() << "+";
          d1->writeOutput(jacobian_output, false, temporary_terms);
          jacobian_output << ";" << endl;
        }
    }

  // Write Hessian w.r. to endogenous only
  if (computeStaticHessian)
    for(second_derivatives_type::const_iterator it = second_derivatives.begin();
        it != second_derivatives.end(); it++)
      {
        int eq = it->first.first;
        int var1 = it->first.second.first;
        int var2 = it->first.second.second;
        NodeID d2 = it->second;

        // Keep only derivatives w.r. to endogenous variables
        if (variable_table.getType(var1) == eEndogenous
            && variable_table.getType(var2) == eEndogenous)
          {
            int id1 = variable_table.getSymbolID(var1);
            int id2 = variable_table.getSymbolID(var2);

            int col_nb = id1*symbol_table.endo_nbr+id2+1;
            int col_nb_sym = id2*symbol_table.endo_nbr+id1+1;

            hessian_output << "  g2" << lpar << eq+1 << ", " << col_nb << rpar << " = ";
            d2->writeOutput(hessian_output, false, temporary_terms);
            hessian_output << ";" << endl;

            // Treating symetric elements
            if (var1 != var2)
              lsymetric <<  "  g2" << lpar << eq+1 << ", " << col_nb_sym << rpar << " = "
                        <<  "g2" << lpar << eq+1 << ", " << col_nb << rpar << ";" << endl;
          }

      }

  // Writing ouputs
  if (offset == 1)
    {
      StaticOutput << "global M_ \n";
      StaticOutput << "if M_.param_nbr > 0\n  params = M_.params;\nend\n";

      StaticOutput << "  residual = zeros( " << equations.size() << ", 1);\n";
      StaticOutput << "\n\t"+interfaces::comment()+"\n\t"+interfaces::comment();
      StaticOutput << "Model equations\n\t";
      StaticOutput << interfaces::comment() + "\n\n";
      StaticOutput << model_output.str();
      StaticOutput << "if ~isreal(residual)\n";
      StaticOutput << "  residual = real(residual)+imag(residual).^2;\n";
      StaticOutput << "end\n";
      StaticOutput << "if nargout >= 2,\n";
      StaticOutput << "  g1 = " <<
        "zeros(" << equations.size() << ", " <<
        symbol_table.endo_nbr << ");\n" ;
      StaticOutput << "\n\t"+interfaces::comment()+"\n\t"+interfaces::comment();
      StaticOutput << "Jacobian matrix\n\t";
      StaticOutput << interfaces::comment() + "\n\n";
      StaticOutput << jacobian_output.str();
      StaticOutput << "  if ~isreal(g1)\n";
      StaticOutput << "    g1 = real(g1)+2*imag(g1);\n";
      StaticOutput << "  end\n";
      StaticOutput << "end\n";
      if (computeStaticHessian)
        {
          StaticOutput << "if nargout >= 3,\n";
          // Writing initialization instruction for matrix g2
          int ncols = symbol_table.endo_nbr * symbol_table.endo_nbr;
          StaticOutput << "  g2 = " <<
            "sparse([],[],[]," << equations.size() << ", " << ncols << ", " <<
            5*ncols << ");\n";
          StaticOutput << "\n\t"+interfaces::comment()+"\n\t"+interfaces::comment();
          StaticOutput << "Hessian matrix\n\t";
          StaticOutput << interfaces::comment() + "\n\n";
          StaticOutput << hessian_output.str() << lsymetric.str();
          StaticOutput << "end;\n";
        }
    }
  else
    {
      StaticOutput << "void Static(double *y, double *x, double *residual, double *g1)\n";
      StaticOutput << "{\n";
      StaticOutput << "  double lhs, rhs;\n\n";
      // Writing residual equations
      StaticOutput << "  /* Residual equations */\n";
      StaticOutput << "  if (residual == NULL) return;\n";
      StaticOutput << " {\n";
      StaticOutput << model_output.str();
      // Writing Jacobian
      StaticOutput << "   /* Jacobian for endogenous variables without lag */\n";
      StaticOutput << "   if (g1 == NULL) return;\n";
      StaticOutput << " {\n";
      StaticOutput << jacobian_output.str();
      StaticOutput << "  }\n";
      StaticOutput << " }\n";
      StaticOutput << "}\n\n";
    }
}

void
ModelTree::writeDynamicModel(ostream &DynamicOutput) const
{
  ostringstream lsymetric;       // Used when writing symetric elements in Hessian
  ostringstream model_output;    // Used for storing model equations
  ostringstream jacobian_output; // Used for storing jacobian equations
  ostringstream hessian_output;  // Used for storing Hessian equations
  ostringstream third_derivatives_output;

  writeTemporaryTerms(model_output, true);

  writeLocalParameters(model_output, true);

  writeModelEquations(model_output, true);

  int nrows = equations.size();
  int nvars;
  if (computeJacobianExo)
    nvars = variable_table.get_dyn_var_nbr();
  else
    nvars = variable_table.var_endo_nbr;
  int nvars_sq = nvars * nvars;

  // Writing Jacobian
  if (computeJacobian || computeJacobianExo)
    for(first_derivatives_type::const_iterator it = first_derivatives.begin();
        it != first_derivatives.end(); it++)
      {
        int eq = it->first.first;
        int var = it->first.second;
        NodeID d1 = it->second;

        if (computeJacobianExo || variable_table.getType(var) == eEndogenous)
          {
            ostringstream g1;
            g1 << "  g1" << lpar << eq + 1 << ", " << variable_table.getSortID(var) + 1 << rpar;

            jacobian_output << g1.str() << "=" << g1.str() << "+";
            d1->writeOutput(jacobian_output, true, temporary_terms);
            jacobian_output << ";" << endl;
          }
      }

  // Writing Hessian
  if (computeHessian)
    for(second_derivatives_type::const_iterator it = second_derivatives.begin();
        it != second_derivatives.end(); it++)
      {
        int eq = it->first.first;
        int var1 = it->first.second.first;
        int var2 = it->first.second.second;
        NodeID d2 = it->second;

        int id1 = variable_table.getSortID(var1);
        int id2 = variable_table.getSortID(var2);

        int col_nb = id1*nvars+id2+1;
        int col_nb_sym = id2*nvars+id1+1;

        hessian_output << "  g2" << lpar << eq+1 << ", " << col_nb << rpar << " = ";
        d2->writeOutput(hessian_output, true, temporary_terms);
        hessian_output << ";" << endl;

        // Treating symetric elements
        if (id1 != id2)
          lsymetric <<  "  g2" << lpar << eq+1 << ", " << col_nb_sym << rpar << " = "
                    <<  "g2" << lpar << eq+1 << ", " << col_nb << rpar << ";" << endl;
      }

  // Writing third derivatives
  if (computeThirdDerivatives)
    for(third_derivatives_type::const_iterator it = third_derivatives.begin();
        it != third_derivatives.end(); it++)
      {
        int eq = it->first.first;
        int var1 = it->first.second.first;
        int var2 = it->first.second.second.first;
        int var3 = it->first.second.second.second;
        NodeID d3 = it->second;

        int id1 = variable_table.getSortID(var1);
        int id2 = variable_table.getSortID(var2);
        int id3 = variable_table.getSortID(var3);

        // Reference column number for the g3 matrix
        int ref_col = id1 * nvars_sq + id2 * nvars + id3 + 1;

        third_derivatives_output << "  g3" << lpar << eq+1 << ", " << ref_col << rpar << " = ";
        d3->writeOutput(third_derivatives_output, true, temporary_terms);
        third_derivatives_output << ";" << endl;

        // Compute the column numbers for the 5 other permutations of (id1,id2,id3) and store them in a set (to avoid duplicates if two indexes are equal)
        set<int> cols;
        cols.insert(id1 * nvars_sq + id3 * nvars + id2 + 1);
        cols.insert(id2 * nvars_sq + id1 * nvars + id3 + 1);
        cols.insert(id2 * nvars_sq + id3 * nvars + id1 + 1);
        cols.insert(id3 * nvars_sq + id1 * nvars + id2 + 1);
        cols.insert(id3 * nvars_sq + id2 * nvars + id1 + 1);

        for(set<int>::iterator it2 = cols.begin(); it2 != cols.end(); it2++)
          if (*it2 != ref_col)
            third_derivatives_output << "  g3" << lpar << eq+1 << ", " << *it2 << rpar << " = "
                                     << "g3" << lpar << eq+1 << ", " << ref_col << rpar
                                     << ";" << endl;
      }

  if (offset == 1)
    {
      DynamicOutput << "global M_ it_\n";
      DynamicOutput << "if M_.param_nbr > 0\n  params =  M_.params;\nend\n";
      DynamicOutput << "\n\t"+interfaces::comment()+"\n\t"+interfaces::comment();
      DynamicOutput << "Model equations\n\t";
      DynamicOutput << interfaces::comment() + "\n\n";
      DynamicOutput << "residual = zeros(" << nrows << ", 1);\n";

      DynamicOutput << model_output.str();

      if (computeJacobian || computeJacobianExo)
        {
          DynamicOutput << "if nargout >= 2,\n";
          // Writing initialization instruction for matrix g1
          DynamicOutput << "  g1 = " <<
            "zeros(" << nrows << ", " << nvars << ");\n" ;
          DynamicOutput << "\n\t"+interfaces::comment()+"\n\t"+interfaces::comment();
          DynamicOutput << "Jacobian matrix\n\t";
          DynamicOutput << interfaces::comment()+"\n\n";
          DynamicOutput << jacobian_output.str();
          DynamicOutput << "end\n";
        }
      if (computeHessian)
        {
          DynamicOutput << "if nargout >= 3,\n";
          // Writing initialization instruction for matrix g2
          int ncols = nvars_sq;
          DynamicOutput << "  g2 = sparse([],[],[]," << nrows << ", " << ncols << ", "
                        << 5*ncols << ");\n";
          DynamicOutput << "\n\t"+interfaces::comment() << "\n\t" << interfaces::comment();
          DynamicOutput << "Hessian matrix\n\t" << interfaces::comment() << "\n\n";
          DynamicOutput << hessian_output.str() << lsymetric.str();
          DynamicOutput << "end;\n";
        }
      if (computeThirdDerivatives)
        {
          DynamicOutput << "if nargout >= 4,\n";
          int ncols = nvars_sq * nvars;
          DynamicOutput << "  g3 = sparse([],[],[]," << nrows << ", " << ncols << ", "
                        << 5*ncols << ");\n";
          DynamicOutput << "\n\t" + interfaces::comment() + "\n\t" + interfaces::comment();
          DynamicOutput << "Third order derivatives\n\t" << interfaces::comment() << "\n\n";
          DynamicOutput << third_derivatives_output.str();
          DynamicOutput << "end;\n";
        }
    }
  else
    {
      DynamicOutput << "void Dynamic(double *y, double *x, double *residual, double *g1, double *g2)\n";
      DynamicOutput << "{\n";
      DynamicOutput << "  double lhs, rhs;\n\n";
      DynamicOutput << "  /* Residual equations */\n";
      DynamicOutput << model_output.str();
      if (computeJacobian || computeJacobianExo)
        {
          DynamicOutput << "  /* Jacobian  */\n";
          DynamicOutput << "  if (g1 == NULL) return;\n";
          DynamicOutput << "  {\n";
          DynamicOutput << jacobian_output.str();
          DynamicOutput << "  }\n";
        }
      if (computeHessian)
        {
          DynamicOutput << "  /* Hessian for endogenous and exogenous variables */\n";
          DynamicOutput << "  if (g2 == NULL) return;\n";
          DynamicOutput << "   {\n";
          DynamicOutput << hessian_output.str() << lsymetric.str();
          DynamicOutput << "   }\n";
        }
      DynamicOutput << "}\n\n";
    }
}

void
ModelTree::writeOutput(ostream &output) const
{
  /* Writing initialisation for M_.lead_lag_incidence matrix
     M_.lead_lag_incidence is a matrix with as many columns as there are
     endogenous variables and as many rows as there are periods in the
     models (nbr of rows = M_.max_lag+M_.max_lead+1)

     The matrix elements are equal to zero if a variable isn't present in the
     model at a given period.
  */
  output << "M_.lead_lag_incidence = [";
  // Loop on endogenous variables
  for (int endoID = 0; endoID < symbol_table.endo_nbr; endoID++)
    {
      output << "\n\t";
      // Loop on periods
      for (int lag = -variable_table.max_endo_lag; lag <= variable_table.max_endo_lead; lag++)
        {
          // Getting name of symbol
          string name = symbol_table.getNameByID(eEndogenous, endoID);
          // and its variableID if exists with current period
          int varID = variable_table.getID(name, lag);
          if (varID >= 0)
            output << " " << variable_table.getPrintIndex(varID) + 1;
          else
            output << " 0";
        }
      output << ";";
    }
  output << "]';\n";

  // Writing initialization for some other variables
  output << "M_.exo_names_orig_ord = [1:" << symbol_table.exo_nbr << "];\n";
  output << "M_.maximum_lag = " << variable_table.max_lag << ";\n";
  output << "M_.maximum_lead = " << variable_table.max_lead << ";\n";
  if (symbol_table.endo_nbr)
    {
      output << "M_.maximum_endo_lag = " << variable_table.max_endo_lag << ";\n";
      output << "M_.maximum_endo_lead = " << variable_table.max_endo_lead << ";\n";
      output << "oo_.steady_state = zeros(" << symbol_table.endo_nbr << ", 1);\n";
    }
  if (symbol_table.exo_nbr)
    {
      output << "M_.maximum_exo_lag = " << variable_table.max_exo_lag << ";\n";
      output << "M_.maximum_exo_lead = " << variable_table.max_exo_lead << ";\n";
      output << "oo_.exo_steady_state = zeros(" << symbol_table.exo_nbr << ", 1);\n";
    }
  if (symbol_table.exo_det_nbr)
    {
      output << "M_.maximum_exo_det_lag = " << variable_table.max_exo_det_lag << ";\n";
      output << "M_.maximum_exo_det_lead = " << variable_table.max_exo_det_lead << ";\n";
      output << "oo_.exo_det_steady_state = zeros(" << symbol_table.exo_det_nbr << ", 1);\n";
    }
  if (symbol_table.recur_nbr)
    {
      output << "M_.maximum_recur_lag = " << variable_table.max_recur_lag << ";\n";
      output << "M_.maximum_recur_lead = " << variable_table.max_recur_lead << ";\n";
      output << "oo_.recur_steady_state = zeros(" << symbol_table.recur_nbr << ", 1);\n";
    }
  if (symbol_table.parameter_nbr)
    output << "M_.params = zeros(" << symbol_table.parameter_nbr << ", 1);\n";
}

void
ModelTree::addEquation(NodeID eq)
{
  BinaryOpNode *beq = dynamic_cast<BinaryOpNode *>(eq);

  if (beq == NULL || beq->op_code != oEqual)
    {
      cerr << "ModelTree::addEquation: you didn't provide an equal node!" << endl;
      exit(-1);
    }

  equations.push_back(beq);
}

void
ModelTree::checkPass() const
{
  // Exit if there is no equation in model file
  if (equations.size() == 0)
    {
      cerr << "No equation found in model file" << endl;
      exit(-1);
    }
}

void
ModelTree::computingPass()
{
  cout << equations.size() << " equation(s) found" << endl;

  // Sorting variable table
  variable_table.Sort();

  if (offset == 1)
    {
      min_cost = 40 * 90;
      lpar = '(';
      rpar = ')';
    }
  else
    {
      min_cost = 40 * 4;
      lpar = '[';
      rpar = ']';
    }

  // Determine derivation order
  int order = 1;
  if (computeThirdDerivatives)
    order = 3;
  else if (computeHessian || computeStaticHessian)
    order = 2;

  // Launch computations
  derive(order);
  computeTemporaryTerms(order);
}

void
ModelTree::writeStaticFile(const string &basename) const
{
  if (offset)
    writeStaticMFile(basename + "_static");
  else
    writeStaticCFile(basename + "_static");
}

void
ModelTree::writeDynamicFile(const string &basename) const
{
  if (offset)
    writeDynamicMFile(basename + "_dynamic");
  else
    writeDynamicCFile(basename + "_dynamic");
}
