#include <iostream>
#include <fstream>
#include <sstream>

#include <math.h>

#include "ModelTree.hh"
#include "Interface.hh"

#include "Model_Graph.hh"
#include "SymbolGaussElim.hh"

ModelTree::ModelTree(SymbolTable &symbol_table_arg,
                     NumericalConstants &num_constants_arg) :
  DataTree(symbol_table_arg, num_constants_arg),
  computeJacobian(false),
  computeJacobianExo(false),
  computeHessian(false),
  computeStaticHessian(false),
  computeThirdDerivatives(false),
  block_triangular(symbol_table_arg)
{
}

int
ModelTree::equation_number() const
{
  return(equations.size());
}

void
ModelTree::GetDerivatives(ostream &output, int eq, int var, int lag, bool is_dynamic,
                          const temporary_terms_type &temporary_terms) const
{
  first_derivatives_type::const_iterator it=first_derivatives.find(make_pair(eq,variable_table.getmVariableSelector(var,lag)));
  if (it != first_derivatives.end())
    (it->second)->writeOutput(output, is_dynamic, temporary_terms, offset);
  else
    output << 0;
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

  bool offs = false;

  if (temporary_terms.size() > 0 && offset == 2)
    {
      output << "double\n";
      offs = true;
    }

  for(temporary_terms_type::const_iterator it = temporary_terms.begin();
      it != temporary_terms.end(); it++)
    {
      (*it)->writeOutput(output, is_dynamic, temporary_terms, offset);
      output << " = ";

      (*it)->writeOutput(output, is_dynamic, tt2, offset);

      // Insert current node into tt2
      tt2.insert(*it);

      if (offs)
        output << ",\n";
      else
        output << ";" << endl;
    }
  if (offs)
    output << ";\n";
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
      value->writeOutput(output, is_dynamic, temporary_terms, offset);
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
      lhs->writeOutput(output, is_dynamic, temporary_terms, offset);
      output << ";" << endl;

      NodeID rhs = eq_node->arg2;
      output << "rhs =";
      rhs->writeOutput(output, is_dynamic, temporary_terms, offset);
      output << ";" << endl;

      output << "residual" << lpar << eq + 1 << rpar << "= lhs-rhs;" << endl;
    }
}

void
ModelTree::computeTemporaryTermsOrdered(int order, Model_Block *ModelBlock)
{
  map<NodeID, int> reference_count, first_occurence;
  int i, j, m, eq, var, lag/*, prev_size=0*/;
  temporary_terms_type vect;
  ostringstream tmp_output;
  BinaryOpNode *eq_node;
  NodeID lhs, rhs;
  first_derivatives_type::const_iterator it;
  ostringstream tmp_s;

  temporary_terms.clear();
  for(j = 0;j < ModelBlock->Size;j++)
    {
      if (ModelBlock->Block_List[j].Size==1)
        {
          eq_node = equations[ModelBlock->Block_List[j].Equation[0]];
          lhs = eq_node->arg1;
          rhs = eq_node->arg2;
          tmp_output.str("");
          lhs->writeOutput(tmp_output, true, temporary_terms, offset);
          tmp_s << "y[Per_y_+" << ModelBlock->Block_List[j].Variable[0] << "]";
          if (tmp_output.str()==tmp_s.str())
            {
              if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_SIMPLE)
                ModelBlock->Block_List[j].Simulation_Type=EVALUATE_BACKWARD;
              else if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_FOREWARD_SIMPLE)
                ModelBlock->Block_List[j].Simulation_Type=EVALUATE_FOREWARD;
            }
        }
      for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
        {
          eq_node = equations[ModelBlock->Block_List[j].Equation[i]];
          eq_node->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock);
        }
      if (ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD
          && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FOREWARD)
        {
          if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE ||
              ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_SIMPLE)
            {
              for(m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
                {
                  lag=m-ModelBlock->Block_List[j].Max_Lag;
                  for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                    {
                      eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                      var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                      it=first_derivatives.find(make_pair(eq,variable_table.getmVariableSelector(var,lag)));
                      it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock);
                    }
                }
            }
          else if (ModelBlock->Block_List[j].Simulation_Type!=SOLVE_BACKWARD_SIMPLE
                   && ModelBlock->Block_List[j].Simulation_Type!=SOLVE_FOREWARD_SIMPLE)
            {
              m=ModelBlock->Block_List[j].Max_Lag;
              for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                {
                  eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                  var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                  it=first_derivatives.find(make_pair(eq,variable_table.getmVariableSelector(var,0)));
                  it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock);
                }
            }
          else
            {
              eq=ModelBlock->Block_List[j].Equation[0];
              var=ModelBlock->Block_List[j].Variable[0];
              it=first_derivatives.find(make_pair(eq,variable_table.getmVariableSelector(var,0)));
              it->second->computeTemporaryTerms(reference_count, temporary_terms, first_occurence, j, ModelBlock);
            }
        }
    }
  if (order == 2)
    for(second_derivatives_type::iterator it = second_derivatives.begin();
        it != second_derivatives.end(); it++)
      it->second->computeTemporaryTerms(reference_count, temporary_terms);
}


void
ModelTree::writeModelEquationsOrdered(ostream &output, bool is_dynamic, Model_Block *ModelBlock) const
{
  int i,j,k,m;
  string sModel, tmp_s;
  ostringstream tmp_output;
  NodeID lhs, rhs;
  BinaryOpNode *eq_node;
  bool OK, lhs_rhs_done, skip_the_head;
  ostringstream Uf[symbol_table.endo_nbr];
  map<NodeID, int> reference_count;
  int prev_Simulation_Type=-1;
  temporary_terms_type::const_iterator it_temp=temporary_terms.begin();
  //----------------------------------------------------------------------
  //Temporary variables décalaration
  OK=true;
  for(temporary_terms_type::const_iterator it = temporary_terms.begin();
      it != temporary_terms.end(); it++)
    {
      if (OK)
        OK=false;
      else
        tmp_output << ", ";

      (*it)->writeOutput(tmp_output, is_dynamic, temporary_terms, 0);

      tmp_output << "[" << block_triangular.periods + variable_table.max_lag+variable_table.max_lead << "]";
    }
  if (tmp_output.str().length()>0)
    {
      output << "double " << tmp_output.str() << ";\n\n";
    }
  //For each block
  for(j = 0;j < ModelBlock->Size;j++)
    {
      //For a block composed of a single equation determines wether we have to evaluate or to solve the equation
      if (ModelBlock->Block_List[j].Size==1)
        {
          lhs_rhs_done=true;
          eq_node = equations[ModelBlock->Block_List[j].Equation[0]];
          lhs = eq_node->arg1;
          rhs = eq_node->arg2;
          tmp_output.str("");
          lhs->writeOutput(tmp_output, is_dynamic, temporary_terms, offset);
        }
      else
        lhs_rhs_done=false;
      if (prev_Simulation_Type==ModelBlock->Block_List[j].Simulation_Type
          && (ModelBlock->Block_List[j].Simulation_Type==EVALUATE_BACKWARD
              ||ModelBlock->Block_List[j].Simulation_Type==EVALUATE_FOREWARD ))
        skip_the_head=true;
      else
        skip_the_head=false;
      if (!skip_the_head)
        {
          if (j>0)
            output << "}\n\n";
          output << "void Dynamic" << j+1 << "(double *y, double *x, double *residual, double *g1, double *g2)\n";
          output << "{\n";
          output << "  ////////////////////////////////////////////////////////////////////////\n" <<
            "  //" << string("                     Block ").substr(int(log10(j + 1))) << j + 1 << " " << BlockTriangular::BlockType0(ModelBlock->Block_List[j].Type) <<
            "          //\n" <<
            "  //                     Simulation type ";
          output << BlockTriangular::BlockSim(ModelBlock->Block_List[j].Simulation_Type) << "  //\n" <<
            "  ////////////////////////////////////////////////////////////////////////\n";
        }
      //The Temporary terms
      temporary_terms_type tt2;
      if (ModelBlock->Block_List[j].Temporary_terms->size())
        output << "  //Temporary variables\n";
      i=0;
      for(temporary_terms_type::const_iterator it = ModelBlock->Block_List[j].Temporary_terms->begin();
          it != ModelBlock->Block_List[j].Temporary_terms->end(); it++)
        {
          output << "  ";
          (*it)->writeOutput(output, is_dynamic, temporary_terms, offset);
          output << " = ";
          (*it)->writeOutput(output, is_dynamic, tt2, offset);
          // Insert current node into tt2
          tt2.insert(*it);
          output << ";" << endl;
          i++;
        }
      // The equations
      for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
        {
          sModel = symbol_table.getNameByID(eEndogenous, ModelBlock->Block_List[j].Variable[i]) ;
          ModelBlock->Block_List[j].Variable_Sorted[i] = variable_table.getID(sModel, 0);
          output << "  //equation " << ModelBlock->Block_List[j].Equation[i] << " variable : " <<
            sModel << " (" << ModelBlock->Block_List[j].Variable[i] << ")\n";
          if (!lhs_rhs_done)
            {
              eq_node = equations[ModelBlock->Block_List[j].Equation[i]];
              lhs = eq_node->arg1;
              rhs = eq_node->arg2;
              tmp_output.str("");
              lhs->writeOutput(tmp_output, is_dynamic, temporary_terms, offset);
            }
          output << "  ";
          switch(ModelBlock->Block_List[j].Simulation_Type)
            {
            case EVALUATE_BACKWARD:
            case EVALUATE_FOREWARD:
              output << tmp_output.str();
              output << " = ";
              rhs->writeOutput(output, is_dynamic, temporary_terms, offset);
              output << ";\n";
              break;
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FOREWARD_COMPLETE:
              Uf[ModelBlock->Block_List[j].Equation[i]] << "  u[" << i << "] = residual[" << i << "]";
              goto end;
            case SOLVE_TWO_BOUNDARIES_COMPLETE:
              Uf[ModelBlock->Block_List[j].Equation[i]] << "  u[" << i << "+Per_u_] = residual[" << i << "]";
              goto end;
            default:
            end:
              output << "residual[" << i << "] = (";
              output << tmp_output.str();
              output << ") - (";
              rhs->writeOutput(output, is_dynamic, temporary_terms, offset);
              output << ");\n";
            }
        }
      // The Jacobian if we have to solve the block
      if (ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_BACKWARD
          && ModelBlock->Block_List[j].Simulation_Type!=EVALUATE_FOREWARD)
        {
          output << "  /* Jacobian  */\n";
          switch(ModelBlock->Block_List[j].Simulation_Type)
            {
            case SOLVE_BACKWARD_SIMPLE:
            case SOLVE_FOREWARD_SIMPLE:
              output << "  g1[0]=";
              GetDerivatives(output, ModelBlock->Block_List[j].Equation[0], ModelBlock->Block_List[j].Variable[0], 0, true, temporary_terms);
              output << "; /* variable=" <<  symbol_table.getNameByID(eEndogenous, ModelBlock->Block_List[j].Variable[0])
                     <<"(" << variable_table.getLag(variable_table.getSymbolID(ModelBlock->Block_List[j].Variable[0])) << ") " << ModelBlock->Block_List[j].Variable[0]
                     << ", equation=" <<  ModelBlock->Block_List[j].Equation[0] << "*/\n";
              break;
            case SOLVE_BACKWARD_COMPLETE:
            case SOLVE_FOREWARD_COMPLETE:
              m=ModelBlock->Block_List[j].Max_Lag;
              for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                {
                  int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                  int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                  int u=ModelBlock->Block_List[j].IM_lead_lag[m].u[i];
                  int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                  int varr=ModelBlock->Block_List[j].IM_lead_lag[m].Var[i];
                  Uf[ModelBlock->Block_List[j].Equation[eqr]] << "-u[" << u << "]*y[Per_y_+" << var << "]";
                  output << "  u[" << u << "] = g1[" << eqr << "*" << ModelBlock->Block_List[j].Size << "+" << varr << "] = ";
                  GetDerivatives(output, eq, var, 0,true, temporary_terms);
                  output << "; // variable=" <<  symbol_table.getNameByID(eEndogenous, var)
                         <<"(" << variable_table.getLag(variable_table.getSymbolID(var))<< ") " << var
                         << ", equation=" <<  eq << "\n";
                }
              for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
                output << Uf[ModelBlock->Block_List[j].Equation[i]].str() << ";\n";
              break;
            case SOLVE_TWO_BOUNDARIES_COMPLETE:
              for(m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
                {
                  k=m-ModelBlock->Block_List[j].Max_Lag;
                  for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                    {
                      int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                      int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                      int u=ModelBlock->Block_List[j].IM_lead_lag[m].u[i];
                      int eqr=ModelBlock->Block_List[j].IM_lead_lag[m].Equ[i];
                      if (k==0)
                        Uf[ModelBlock->Block_List[j].Equation[eqr]] << "-u[" << u << "+Per_u_]*y[Per_y_+" << var << "]";
                      else if (k>0)
                        Uf[ModelBlock->Block_List[j].Equation[eqr]] << "-u[" << u << "+Per_u_]*y[(it_+" << k << ")*y_size+" << var << "]";
                      else if (k<0)
                        Uf[ModelBlock->Block_List[j].Equation[eqr]] << "-u[" << u << "+Per_u_]*y[(it_" << k << ")*y_size+" << var << "]";
                      output << "  u[" << u << "+Per_u_] = ";
                      GetDerivatives(output, eq, var, k,true,temporary_terms);
                      output << "; // variable=" <<  symbol_table.getNameByID(eEndogenous, var)
                             <<"(" << k << ") " << var
                             << ", equation=" <<  eq << "\n";
                    }
                }
              for(i = 0;i < ModelBlock->Block_List[j].Size;i++)
                output << Uf[ModelBlock->Block_List[j].Equation[i]].str() << ";\n";
              break;
            }
        }
      prev_Simulation_Type=ModelBlock->Block_List[j].Simulation_Type;
    }
  output << "}\n\n";
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

  writeDynamicModel(mDynamicModelFile, block_triangular.ModelBlock);

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
  mStaticModelFile << "     mexPrintf(\"Global variable not found : \");\n";
  mStaticModelFile << "     mexErrMsgTxt(\"M_ \\n\");\n";
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
  string filename;
  ofstream mDynamicModelFile;
  string tmp_s;
  int i, j;
  if (offset == 2)
    {
      if (compiler==GCC_COMPILE)
        filename = dynamic_basename + ".hh";
      else
        filename = dynamic_basename + ".h";
      mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
      if (!mDynamicModelFile.is_open())
        {
          cout << "ModelTree::Open : Error : Can't open file " << filename
               << ".h for writing\n";
          exit(-1);
        }
      filename.erase(filename.end() - 2, filename.end());
      tmp_s = filename;
      j = tmp_s.size();
      for(i = 0;i < j;i++)
        if ((tmp_s[i] == '\\') || (tmp_s[i] == '.') || (tmp_s[i] == ':'))
          tmp_s[i] = '_';
      mDynamicModelFile << "#ifndef " << tmp_s << "\n";
      mDynamicModelFile << "#define " << tmp_s << "\n";
      if (compiler==GCC_COMPILE)
        {
          mDynamicModelFile << "typedef struct IM_compact\n";
          mDynamicModelFile << "{\n";
          mDynamicModelFile << "  int size, u_init, u_finish, nb_endo;\n";
          mDynamicModelFile << "  int *u, *Var, *Equ, *Var_Index, *Equ_Index, *Var_dyn_Index;\n";
          mDynamicModelFile << "};\n";
          mDynamicModelFile << "typedef struct Variable_l\n";
          mDynamicModelFile << "{\n";
          mDynamicModelFile << "  int* Index;\n";
          mDynamicModelFile << "};\n";
          mDynamicModelFile << "typedef struct tBlock\n";
          mDynamicModelFile << "{\n";
          mDynamicModelFile << "    int Size, Sized, Type, Max_Lead, Max_Lag, Simulation_Type, /*icc1_size,*/ Nb_Lead_Lag_Endo;\n";
          mDynamicModelFile << "    int *Variable, *dVariable, *Equation/*, *icc1, *ics*/;\n";
          mDynamicModelFile << "    int *variable_dyn_index, *variable_dyn_leadlag;\n";
          mDynamicModelFile << "    IM_compact *IM_lead_lag;\n";
          mDynamicModelFile << "};\n";
          mDynamicModelFile << "\n";
          mDynamicModelFile << "typedef struct tModel_Block\n";
          mDynamicModelFile << "{\n";
          mDynamicModelFile << "    int Size;\n";
          mDynamicModelFile << "    tBlock * List;\n";
          mDynamicModelFile << "};\n";
          mDynamicModelFile << "\n";
        }
      else
        {
          mDynamicModelFile << "typedef struct IM_compact\n";
          mDynamicModelFile << "{\n";
          mDynamicModelFile << "  int size, u_init, u_finish, nb_endo;\n";
          mDynamicModelFile << "  int *u, *Var, *Equ, *Var_Index, *Equ_Index, *Var_dyn_Index;\n";
          mDynamicModelFile << "} IM_compact;\n";
          mDynamicModelFile << "typedef struct Variable_l\n";
          mDynamicModelFile << "{\n";
          mDynamicModelFile << "  int* Index;\n";
          mDynamicModelFile << "} Variable_l;\n";
          mDynamicModelFile << "typedef struct tBlock\n";
          mDynamicModelFile << "{\n";
          mDynamicModelFile << "    int Size, Sized, Type, Max_Lead, Max_Lag, Simulation_Type, /*icc1_size,*/ Nb_Lead_Lag_Endo;\n";
          mDynamicModelFile << "    int *Variable, *dVariable, *Equation/*, *icc1, *ics*/;\n";
          mDynamicModelFile << "    int *variable_dyn_index, *variable_dyn_leadlag;\n";
          mDynamicModelFile << "    IM_compact *IM_lead_lag;\n";
          mDynamicModelFile << "} tBlock;\n";
          mDynamicModelFile << "\n";
          mDynamicModelFile << "typedef struct tModel_Block\n";
          mDynamicModelFile << "{\n";
          mDynamicModelFile << "    int Size;\n";
          mDynamicModelFile << "    tBlock * List;\n";
          mDynamicModelFile << "} tModel_Block;\n";
          mDynamicModelFile << "\n";
          mDynamicModelFile << "double *u, slowc, max_res, res2, res1;\n";
          mDynamicModelFile << "double *params;\n";
          mDynamicModelFile << "int it_,Per_u_;\n";
          mDynamicModelFile << "bool cvg;\n";
          mDynamicModelFile << "int nb_row_x;\n";
          mDynamicModelFile << "int y_kmin, y_kmax,periods, x_size, y_size, u_size, maxit_;\n";
          mDynamicModelFile << "double *y=NULL, *x=NULL, *r=NULL, *g1=NULL, *g2=NULL, solve_tolf, dynaretol;\n";
          mDynamicModelFile << "pctimer_t t0, t1;\n";
        }
      mDynamicModelFile << "const int UNKNOWN=" << UNKNOWN << ";\n";
      mDynamicModelFile << "const int EVALUATE_FOREWARD=" << EVALUATE_FOREWARD << ";\n";
      mDynamicModelFile << "const int EVALUATE_BACKWARD=" << EVALUATE_BACKWARD << ";\n";
      mDynamicModelFile << "const int SOLVE_FOREWARD_SIMPLE=" << SOLVE_FOREWARD_SIMPLE << ";\n";
      mDynamicModelFile << "const int SOLVE_BACKWARD_SIMPLE=" << SOLVE_BACKWARD_SIMPLE << ";\n";
      mDynamicModelFile << "const int SOLVE_TWO_BOUNDARIES_SIMPLE=" << SOLVE_TWO_BOUNDARIES_SIMPLE << ";\n";
      mDynamicModelFile << "const int SOLVE_FOREWARD_COMPLETE=" << SOLVE_FOREWARD_COMPLETE << ";\n";
      mDynamicModelFile << "const int SOLVE_BACKWARD_COMPLETE=" << SOLVE_BACKWARD_COMPLETE << ";\n";
      mDynamicModelFile << "const int SOLVE_TWO_BOUNDARIES_COMPLETE=" << SOLVE_TWO_BOUNDARIES_COMPLETE << ";\n";
      mDynamicModelFile << "#endif\n";
      mDynamicModelFile.close();
    }
  if (offset == 1||(offset==2 && compiler==LCC_COMPILE))
    filename = dynamic_basename + ".c";
  else
    filename = dynamic_basename + ".cc";
  mDynamicModelFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDynamicModelFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(-1);
    }
  mDynamicModelFile << "/*\n";
  mDynamicModelFile << " * " << filename << " : Computes dynamic model for Dynare\n";
  mDynamicModelFile << " *\n";
  mDynamicModelFile << " * Warning : this file is generated automatically by Dynare\n";
  mDynamicModelFile << " *           from model file (.mod)\n\n";
  mDynamicModelFile << " */\n";
  if (offset == 2)
    {
      if (compiler==GCC_COMPILE)
        {
          mDynamicModelFile << "#include \"" << dynamic_basename.c_str() << ".hh\"\n";
          mDynamicModelFile << "#include \"simulate.cc\"\n";
        }
      else
        {
          mDynamicModelFile << "#include <math.h>\n";
          mDynamicModelFile << "#include <stdio.h>\n";
          mDynamicModelFile << "#include <string.h>\n";
          mDynamicModelFile << "#include \"pctimer_h.h\"\n";
          mDynamicModelFile << "#include \"mex.h\" /* The Last include file*/\n";
          mDynamicModelFile << "#include \"" << dynamic_basename.c_str() << ".h\"\n";
          mDynamicModelFile << "#include \"simulate.h\"\n";
        }
    }
  if (offset == 2)
    mDynamicModelFile << "//#define DEBUG\n";
  if (offset == 0)
    {
      // A flobal variable for model parameters
      mDynamicModelFile << "double *params;\n";
      // A global variable for it_
      mDynamicModelFile << "int it_;\n";
      mDynamicModelFile << "int nb_row_x;\n";
    }
  else
    {
      if (compiler==GCC_COMPILE)
        {
          mDynamicModelFile << "\n";
        }
    }
  // Writing the function body
  if (offset == 2)
    {
      writeDynamicModel(mDynamicModelFile, block_triangular.ModelBlock);
      SaveCFiles(block_triangular.ModelBlock, block_triangular.file_name, mDynamicModelFile);
    }
  else
    writeDynamicModel(mDynamicModelFile,block_triangular.ModelBlock);

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
  mDynamicModelFile << "      mexPrintf(\"Global variable not found : \");\n";
  mDynamicModelFile << "      mexErrMsgTxt(\"M_ \\n\");\n";
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
          d1->writeOutput(jacobian_output, false, temporary_terms, offset);
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
            d2->writeOutput(hessian_output, false, temporary_terms, offset);
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

string
ModelTree::reform(const string name1) const
{
  string name=name1;
  int pos = name.find("\\", 0);
  while(pos >= 0)
    {
      if (name.substr(pos + 1, 1) != "\\")
        {
          name = name.insert(pos, "\\");
          pos++;
        }
      pos++;
      pos = name.find("\\", pos);
    }
  return (name);
}

void
ModelTree::SaveCFiles(Model_Block* ModelBlock, string Model_file_name, ofstream &mDynamicModelFile) const
{
  int i, j, k, Nb_SGE=0;
  bool printed = false, skip_head, open_par=false;
  SymbolicGaussElimination SGE;

  if (mDynamicModelFile.is_open() && (computeJacobian || computeJacobianExo || computeHessian))
    {
      if (offset == 2)
        {
          mDynamicModelFile << "void Dynamic_Init(tModel_Block *Model_Block)\n";
          mDynamicModelFile << " {\n";
          mDynamicModelFile << "   int i;\n";
          int prev_Simulation_Type=-1;
          for(i = 0;i < block_triangular.ModelBlock->Size;i++)
            {
              k = block_triangular.ModelBlock->Block_List[i].Simulation_Type;
              if (prev_Simulation_Type==k &&
                  (k==EVALUATE_FOREWARD || k==EVALUATE_BACKWARD))
                skip_head=true;
              else
                skip_head=false;
              if ((k == EVALUATE_FOREWARD) && (block_triangular.ModelBlock->Block_List[i].Size))
                {
                  if (!skip_head)
                    {
                      if (open_par)
                        {
                          mDynamicModelFile << "#endif\n";
                          mDynamicModelFile << "      }\n";
                        }
                      mDynamicModelFile << "    for(it_=y_kmin;it_<periods+y_kmin;it_++)\n";
                      mDynamicModelFile << "      {\n";
                      mDynamicModelFile << "        Per_y_=it_*y_size;\n";
                      mDynamicModelFile << "        Dynamic" << i + 1 << "(y, x, r, g1, g2);\n";
                      mDynamicModelFile << "#ifdef DEBUG\n";
                    }
                  for(j = 0;j < block_triangular.ModelBlock->Block_List[i].Size;j++)
                    mDynamicModelFile << "        mexPrintf(\"y[%d, %d]=%f \\n\",it_," << block_triangular.ModelBlock->Block_List[i].Variable[j] << ",y[it_," << block_triangular.ModelBlock->Block_List[i].Variable[j] << "]);\n";
                  open_par=true;
                }
              else if ((k == EVALUATE_BACKWARD) && (block_triangular.ModelBlock->Block_List[i].Size))
                {
                  if (!skip_head)
                    {
                      if (open_par)
                        {
                          mDynamicModelFile << "#endif\n";
                          mDynamicModelFile << "      }\n";
                        }
                      mDynamicModelFile << "    for(it_=periods+y_kmin;it_>y_kmin;it_--)\n";
                      mDynamicModelFile << "      {\n";
                      mDynamicModelFile << "        Per_y_=it_*y_size;\n";
                      mDynamicModelFile << "        Dynamic" << i + 1 << "(y, x, r, g1, g2);\n";
                      mDynamicModelFile << "#ifdef DEBUG\n";
                    }
                  for(j = 0;j < block_triangular.ModelBlock->Block_List[i].Size;j++)
                    mDynamicModelFile << "        mexPrintf(\"y[%d, %d]=%f \\n\",it_," << block_triangular.ModelBlock->Block_List[i].Variable[j] << ",y[it_," << block_triangular.ModelBlock->Block_List[i].Variable[j] << "]);\n";
                  open_par=true;
                }
              else if ((k == SOLVE_FOREWARD_SIMPLE) && (block_triangular.ModelBlock->Block_List[i].Size))
                {
                  if (open_par)
                    {
                      mDynamicModelFile << "#endif\n";
                      mDynamicModelFile << "      }\n";
                    }
                  open_par=false;
                  mDynamicModelFile << "    g1=(double*)mxMalloc(" << block_triangular.ModelBlock->Block_List[i].Size*block_triangular.ModelBlock->Block_List[i].Size << "*sizeof(double));\n";
                  mDynamicModelFile << "    r=(double*)mxMalloc(" << block_triangular.ModelBlock->Block_List[i].Size << "*sizeof(double));\n";
                  mDynamicModelFile << "    for(it_=y_kmin;it_<periods+y_kmin;it_++)\n";
                  mDynamicModelFile << "      {\n";
                  mDynamicModelFile << "        cvg=false;\n";
                  mDynamicModelFile << "        iter=0;\n";
                  mDynamicModelFile << "        Per_y_=it_*y_size;\n";
                  mDynamicModelFile << "        while(!((cvg)||(iter>maxit_)))\n";
                  mDynamicModelFile << "          {\n";
                  mDynamicModelFile << "            Dynamic" << i + 1 << "(y, x, r, g1, g2);\n";
                  mDynamicModelFile << "            y[Per_y_+" << block_triangular.ModelBlock->Block_List[i].Variable[0] << "] += -r[0]/g1[0];\n";
                  mDynamicModelFile << "            cvg=((r[0]*r[0])<solve_tolf);\n";
                  mDynamicModelFile << "            iter++;\n";
                  mDynamicModelFile << "          }\n";
                  mDynamicModelFile << "        if (!cvg)\n";
                  mDynamicModelFile << "          {\n";
                  mDynamicModelFile << "            mexPrintf(\"Convergence not achieved in block " << i << ", at time %d after %d iterations\\n\",it_,iter);\n";
                  mDynamicModelFile << "            mexErrMsgTxt(\"End of simulate\");\n";
                  mDynamicModelFile << "          }\n";
                  mDynamicModelFile << "#ifdef DEBUG\n";
                  mDynamicModelFile << "        mexPrintf(\"y[%d, %d]=%f \\n\",it_," << block_triangular.ModelBlock->Block_List[i].Variable[0] << ",y[it_," << block_triangular.ModelBlock->Block_List[i].Variable[0] << "]);\n";
                  mDynamicModelFile << "#endif\n";
                  mDynamicModelFile << "      }\n";
                  mDynamicModelFile << "    mxFree(g1);\n";
                  mDynamicModelFile << "    mxFree(r);\n";
                }
              else if ((k == SOLVE_BACKWARD_SIMPLE) && (block_triangular.ModelBlock->Block_List[i].Size))
                {
                  if (open_par)
                    {
                      mDynamicModelFile << "#endif\n";
                      mDynamicModelFile << "      }\n";
                    }
                  open_par=false;
                  mDynamicModelFile << "    g1=(double*)mxMalloc(" << block_triangular.ModelBlock->Block_List[i].Size*block_triangular.ModelBlock->Block_List[i].Size << "*sizeof(double));\n";
                  mDynamicModelFile << "    r=(double*)mxMalloc(" << block_triangular.ModelBlock->Block_List[i].Size << "*sizeof(double));\n";
                  mDynamicModelFile << "    for(it_=periods+y_kmin;it_>y_kmin;it_--)\n";
                  mDynamicModelFile << "      {\n";
                  mDynamicModelFile << "        cvg=false;\n";
                  mDynamicModelFile << "        iter=0;\n";
                  mDynamicModelFile << "        Per_y_=it_*y_size;\n";
                  mDynamicModelFile << "        while(!((cvg)||(iter>maxit_)))\n";
                  mDynamicModelFile << "          {\n";
                  mDynamicModelFile << "            Dynamic" << i + 1 << "(y, x, r, g1, g2);\n";
                  mDynamicModelFile << "            y[Per_y_+" << block_triangular.ModelBlock->Block_List[i].Variable[0] << "] += -r[0]/g1[0];\n";
                  mDynamicModelFile << "            cvg=((r[0]*r[0])<solve_tolf);\n";
                  mDynamicModelFile << "            iter++;\n";
                  mDynamicModelFile << "          }\n";
                  mDynamicModelFile << "        if (!cvg)\n";
                  mDynamicModelFile << "          {\n";
                  mDynamicModelFile << "            mexPrintf(\"Convergence not achieved in block " << i << ", at time %d after %d iterations\\n\",it_,iter);\n";
                  mDynamicModelFile << "            mexErrMsgTxt(\"End of simulate\");\n";
                  mDynamicModelFile << "          }\n";
                  mDynamicModelFile << "#ifdef DEBUG\n";
                  mDynamicModelFile << "        mexPrintf(\"y[%d, %d]=%f \\n\",it_," << block_triangular.ModelBlock->Block_List[i].Variable[0] << ",y[it_," << block_triangular.ModelBlock->Block_List[i].Variable[0] << "]);\n";
                  mDynamicModelFile << "#endif\n";
                  mDynamicModelFile << "      }\n";
                  mDynamicModelFile << "    mxFree(g1);\n";
                  mDynamicModelFile << "    mxFree(r);\n";
                }
              else if ((k == SOLVE_TWO_BOUNDARIES_SIMPLE) && (block_triangular.ModelBlock->Block_List[i].Size))
                {
                  if (open_par)
                    {
                      mDynamicModelFile << "#endif\n";
                      mDynamicModelFile << "      }\n";
                    }
                  open_par=false;
                  if (!printed)
                    {
                      printed = true;
                    }
                  SGE.SGE_compute(block_triangular.ModelBlock, i, true, Model_file_name, /*mod_param.endo_nbr*/symbol_table.endo_nbr);
                  Nb_SGE++;
#ifdef PRINT_OUT
                  cout << "end of Gaussian elimination\n";
#endif
                  mDynamicModelFile << "    Read_file(\"" << reform(Model_file_name) << "\",periods," <<
                    block_triangular.ModelBlock->Block_List[i].IM_lead_lag[block_triangular.ModelBlock->Block_List[i].Max_Lag + block_triangular.ModelBlock->Block_List[i].Max_Lead].u_finish + 1 << ", " << /*mod_param.endo_nbr*/symbol_table.endo_nbr <<
                    ", " << block_triangular.ModelBlock->Block_List[i].Max_Lag << ", " << block_triangular.ModelBlock->Block_List[i].Max_Lead << ");\n";
                  mDynamicModelFile << "    g1=(double*)mxMalloc(" << block_triangular.ModelBlock->Block_List[i].Size*block_triangular.ModelBlock->Block_List[i].Size << "*sizeof(double));\n";
                  mDynamicModelFile << "    r=(double*)mxMalloc(" << block_triangular.ModelBlock->Block_List[i].Size << "*sizeof(double));\n";
                  if (!block_triangular.ModelBlock->Block_List[i].is_linear)
                    {
                      mDynamicModelFile << "    cvg=false;\n";
                      mDynamicModelFile << "    iter=0;\n";
                      mDynamicModelFile << "    while(!((cvg)||(iter>maxit_)))\n";
                      mDynamicModelFile << "      {\n";
                      mDynamicModelFile << "        res2=0;\n";
                      mDynamicModelFile << "        res1=0;\n";
                      mDynamicModelFile << "        max_res=0;\n";
                      mDynamicModelFile << "        for(it_=y_kmin;it_<periods+y_kmin;it_++)\n";
                      mDynamicModelFile << "          {\n";
                      mDynamicModelFile << "            Per_u_=(it_-y_kmin)*" << block_triangular.ModelBlock->Block_List[i].IM_lead_lag[block_triangular.ModelBlock->Block_List[i].Max_Lag + block_triangular.ModelBlock->Block_List[i].Max_Lead].u_finish + 1 << ";\n";
                      mDynamicModelFile << "            Per_y_=it_*y_size;\n";
                      mDynamicModelFile << "            Dynamic" << i + 1 << "(y, x, r, g1, g2);\n";
                      mDynamicModelFile << "            for(i=0;i<" << block_triangular.ModelBlock->Block_List[i].Size << ";i++)\n";
                      mDynamicModelFile << "              {\n";
                      mDynamicModelFile << "                if (max_res<fabs(r[i]))\n";
                      mDynamicModelFile << "                  max_res=fabs(r[i]);\n";
                      mDynamicModelFile << "                res2+=r[i]*r[i];\n";
                      mDynamicModelFile << "                res1+=fabs(r[i]);\n";
                      mDynamicModelFile << "              }\n";
                      mDynamicModelFile << "          }\n";
                      mDynamicModelFile << "        iter++;\n";
                      mDynamicModelFile << "        cvg=(max_res<solve_tolf);\n";
                      mDynamicModelFile << "        simulate(" << i << ", " << /*mod_param.endo_nbr*/symbol_table.endo_nbr << ", periods, y_kmin, y_kmax);\n";
                      mDynamicModelFile << "      }\n";
                      mDynamicModelFile << "    if (!cvg)\n";
                      mDynamicModelFile << "      {\n";
                      mDynamicModelFile << "        mexPrintf(\"Convergence not achieved in block " << i << ", after %d iterations\\n\",iter);\n";
                      mDynamicModelFile << "        mexErrMsgTxt(\"End of simulate\");\n";
                      mDynamicModelFile << "      }\n";
                    }
                  else
                    {
                      mDynamicModelFile << "    for(it_=y_kmin;it_<periods+y_kmin;it_++)\n";
                      mDynamicModelFile << "      {\n";
                      mDynamicModelFile << "        Per_u_=(it_-y_kmin)*" << block_triangular.ModelBlock->Block_List[i].IM_lead_lag[block_triangular.ModelBlock->Block_List[i].Max_Lag + block_triangular.ModelBlock->Block_List[i].Max_Lead].u_finish + 1 << ";\n";
                      mDynamicModelFile << "        Per_y_=it_*y_size;\n";
                      mDynamicModelFile << "        Dynamic" << i + 1 << "(y, x, r, g1, g2);\n";
                      mDynamicModelFile << "#ifdef PRINT_OUT\n";
                      mDynamicModelFile << "        for(j=0;j<" << block_triangular.ModelBlock->Block_List[i].IM_lead_lag[block_triangular.ModelBlock->Block_List[i].Max_Lag + block_triangular.ModelBlock->Block_List[i].Max_Lead].u_finish + 1 << ";j++)\n";
                      mDynamicModelFile << "          {\n";
                      mDynamicModelFile << "            mexPrintf(\" %f\",u[Per_u_+j]);\n";
                      mDynamicModelFile << "          }\n";
                      mDynamicModelFile << "        mexPrintf(\"\\n\");\n";
                      mDynamicModelFile << "#endif\n";
                      mDynamicModelFile << "      }\n";
                      mDynamicModelFile << "    simulate(" << i << ", " << /*mod_param.endo_nbr*/symbol_table.endo_nbr << ", it_, y_kmin, y_kmax);\n";
                    }
                  mDynamicModelFile << "#ifdef DEBUG\n";
                  mDynamicModelFile << "    for(it_=y_kmin;it_<periods+y_kmin;it_++)\n";
                  mDynamicModelFile << "      {\n";
                  mDynamicModelFile << "        for(i=0;i<Model_Block->List[" << i << "].Size;i++)\n";
                  mDynamicModelFile << "          {";
                  mDynamicModelFile << "            Per_y_=it_*y_size;\n";
                  mDynamicModelFile << "            mexPrintf(\" y[%d, %d]=%f \",it_,Model_Block->List[" << i << "].Variable[i],y[Per_y_+Model_Block->List[" << i << "].Variable[i]]);\n";
                  mDynamicModelFile << "          }";
                  mDynamicModelFile << "        mexPrintf(\" \\n \");\n";
                  mDynamicModelFile << "      }\n";
                  mDynamicModelFile << "#endif\n";
                  mDynamicModelFile << "    mxFree(g1);\n";
                  mDynamicModelFile << "    mxFree(r);\n";
                  mDynamicModelFile << "    mxFree(u);\n";
                  mDynamicModelFile << "    //mexErrMsgTxt(\"Exit from Dynare\");\n";
                }
              else if ((k == SOLVE_FOREWARD_COMPLETE) && (block_triangular.ModelBlock->Block_List[i].Size))
                {
                  if (open_par)
                    {
                      mDynamicModelFile << "#endif\n";
                      mDynamicModelFile << "      }\n";
                    }
                  open_par=false;
                  if (!printed)
                    {
                      printed = true;
                    }
                  SGE.SGE_compute(block_triangular.ModelBlock, i, false, Model_file_name, /*mod_param.endo_nbr*/symbol_table.endo_nbr);
                  Nb_SGE++;
                  mDynamicModelFile << "    Read_file(\"" << reform(Model_file_name) << "\", periods, 0, 0, 0, 0);\n";
                  mDynamicModelFile << "    g1=(double*)mxMalloc(" << block_triangular.ModelBlock->Block_List[i].Size*block_triangular.ModelBlock->Block_List[i].Size << "*sizeof(double));\n";
                  mDynamicModelFile << "    r=(double*)mxMalloc(" << block_triangular.ModelBlock->Block_List[i].Size << "*sizeof(double));\n";
                  mDynamicModelFile << "    for(it_=y_kmin;it_<periods+y_kmin;it_++)\n";
                  mDynamicModelFile << "      {\n";
                  if (!block_triangular.ModelBlock->Block_List[i].is_linear)
                    {
                      mDynamicModelFile << "        cvg=false;\n";
                      mDynamicModelFile << "        iter=0;\n";
                      mDynamicModelFile << "        Per_y_=it_*y_size;\n";
                      mDynamicModelFile << "        while(!((cvg)||(iter>maxit_)))\n";
                      mDynamicModelFile << "          {\n";
                      mDynamicModelFile << "            Dynamic" << i + 1 << "(y, x, r, g1, g2);\n";
                      mDynamicModelFile << "            simulate(" << i << ", " << /*mod_param.endo_nbr*/symbol_table.endo_nbr << ", it_, y_kmin, y_kmax);\n";
                      mDynamicModelFile << "            res2=0;\n";
                      mDynamicModelFile << "            res1=0;\n";
                      mDynamicModelFile << "            max_res=0;\n";
                      mDynamicModelFile << "            for(i=0;i<" << block_triangular.ModelBlock->Block_List[i].Size << ";i++)\n";
                      mDynamicModelFile << "              {\n";
                      mDynamicModelFile << "                if (max_res<fabs(r[i]))\n";
                      mDynamicModelFile << "                  max_res=fabs(r[i]);\n";
                      mDynamicModelFile << "                res2+=r[i]*r[i];\n";
                      mDynamicModelFile << "                res1+=fabs(r[i]);\n";
                      mDynamicModelFile << "              }\n";
                      mDynamicModelFile << "            cvg=(max_res<solve_tolf);\n";
                      mDynamicModelFile << "            iter++;\n";
                      mDynamicModelFile << "          }\n";
                      mDynamicModelFile << "        if (!cvg)\n";
                      mDynamicModelFile << "          {\n";
                      mDynamicModelFile << "            mexPrintf(\"Convergence not achieved in block " << i << ", at time %d after %d iterations\\n\",it_,iter);\n";
                      mDynamicModelFile << "            mexErrMsgTxt(\"End of simulate\");\n";
                      mDynamicModelFile << "          }\n";
                    }
                  else
                    {
                      mDynamicModelFile << "        Dynamic" << i + 1 << "(y, x, r, g1, g2);\n";
                      mDynamicModelFile << "        simulate(" << i << ", " << /*mod_param.endo_nbr*/symbol_table.endo_nbr << ", it_, y_kmin, y_kmax);\n";
                    }
                  mDynamicModelFile << "#ifdef DEBUG\n";
                  mDynamicModelFile << "        for(i=0;i<Model_Block->List[" << i << "].Size;i++)\n";
                  mDynamicModelFile << "          mexPrintf(\" y[%d, %d]=%f \",it_,Model_Block->List[" << i << "].Variable[i],y[it_*" << /*mod_param.endo_nbr*/symbol_table.endo_nbr << "+Model_Block->List[" << i << "].Variable[i]]);\n";
                  mDynamicModelFile << "        mexPrintf(\" \\n \");\n";
                  mDynamicModelFile << "#endif\n";
                  mDynamicModelFile << "      }\n";
                  mDynamicModelFile << "    mxFree(g1);\n";
                  mDynamicModelFile << "    mxFree(r);\n";
                  mDynamicModelFile << "    mxFree(u);\n";
                }
              else if ((k == SOLVE_BACKWARD_COMPLETE) && (block_triangular.ModelBlock->Block_List[i].Size))
                {
                  if (open_par)
                    {
                      mDynamicModelFile << "#endif\n";
                      mDynamicModelFile << "      }\n";
                    }
                  open_par=false;
                  SGE.SGE_compute(block_triangular.ModelBlock, i, false, Model_file_name, /*mod_param.endo_nbr*/symbol_table.endo_nbr);
                  Nb_SGE++;
                  mDynamicModelFile << "    Read_file(\"" << reform(Model_file_name) << "\", periods, 0, 0, 0, 0);\n";
                  mDynamicModelFile << "    g1=(double*)mxMalloc(" << block_triangular.ModelBlock->Block_List[i].Size*block_triangular.ModelBlock->Block_List[i].Size << "*sizeof(double));\n";
                  mDynamicModelFile << "    r=(double*)mxMalloc(" << block_triangular.ModelBlock->Block_List[i].Size << "*sizeof(double));\n";
                  mDynamicModelFile << "    for(it_=periods+y_kmin;it_>y_kmin;it_--)\n";
                  mDynamicModelFile << "      {\n";
                  if (!block_triangular.ModelBlock->Block_List[i].is_linear)
                    {
                      mDynamicModelFile << "        cvg=false;\n";
                      mDynamicModelFile << "        iter=0;\n";
                      mDynamicModelFile << "        Per_y_=it_*y_size;\n";
                      mDynamicModelFile << "        while(!((cvg)||(iter>maxit_)))\n";
                      mDynamicModelFile << "          {\n";
                      mDynamicModelFile << "            Dynamic" << i + 1 << "(y, x, r, g1, g2);\n";
                      mDynamicModelFile << "            simulate(" << i << ", " << /*mod_param.endo_nbr*/symbol_table.endo_nbr << ", it_, y_kmin, y_kmax);\n";
                      mDynamicModelFile << "            res2=0;\n";
                      mDynamicModelFile << "            for(i=0;i<" << block_triangular.ModelBlock->Block_List[i].Size << ";i++)\n";
                      mDynamicModelFile << "              res2+=r[i]*r[i];\n";
                      mDynamicModelFile << "            cvg=(res2<solve_tolf);\n";
                      mDynamicModelFile << "            iter++;\n";
                      mDynamicModelFile << "          }\n";
                      mDynamicModelFile << "        if (!cvg)\n";
                      mDynamicModelFile << "          {\n";
                      mDynamicModelFile << "            mexPrintf(\"Convergence not achieved in block " << i << ", at time %d after %d iterations\\n\",it_,iter);\n";
                      mDynamicModelFile << "            mexErrMsgTxt(\"End of simulate\");\n";
                      mDynamicModelFile << "          }\n";
                    }
                  else
                    {
                      mDynamicModelFile << "        Dynamic" << i + 1 << "(y, x, r, g1, g2);\n";
                      mDynamicModelFile << "        simulate(" << i << ", " << /*mod_param.endo_nbr*/symbol_table.endo_nbr << ", it_, y_kmin, y_kmax);\n";
                    }
                  mDynamicModelFile << "#ifdef DEBUG\n";
                  mDynamicModelFile << "        for(i=0;i<Model_Block->List[" << i << "].Size;i++)\n";
                  mDynamicModelFile << "          mexPrintf(\" y[%d, %d]=%f \",it_,Model_Block->List[" << i << "].Variable[i],y[it_*" << /*mod_param.endo_nbr*/symbol_table.endo_nbr << "+Model_Block->List[" << i << "].Variable[i]]);\n";
                  mDynamicModelFile << "        mexPrintf(\" \\n \");\n";
                  mDynamicModelFile << "#endif\n";
                  mDynamicModelFile << "      }\n";
                  mDynamicModelFile << "    mxFree(g1);\n";
                  mDynamicModelFile << "    mxFree(r);\n";
                  mDynamicModelFile << "    mxFree(u);\n";
                }
              else if ((k == SOLVE_TWO_BOUNDARIES_COMPLETE) && (block_triangular.ModelBlock->Block_List[i].Size))
                {
                  if (open_par)
                    {
                      mDynamicModelFile << "#endif\n";
                      mDynamicModelFile << "      }\n";
                    }
                  open_par=false;
                  if (!printed)
                    {
                      printed = true;
                    }
                  SGE.SGE_compute(block_triangular.ModelBlock, i, true, Model_file_name, /*mod_param.endo_nbr*/symbol_table.endo_nbr);
                  Nb_SGE++;
                  mDynamicModelFile << "    Read_file(\"" << reform(Model_file_name) << "\",periods," <<
                    block_triangular.ModelBlock->Block_List[i].IM_lead_lag[block_triangular.ModelBlock->Block_List[i].Max_Lag + block_triangular.ModelBlock->Block_List[i].Max_Lead].u_finish + 1 << ", " << /*mod_param.endo_nbr*/symbol_table.endo_nbr <<
                    ", " << block_triangular.ModelBlock->Block_List[i].Max_Lag << ", " << block_triangular.ModelBlock->Block_List[i].Max_Lead << ");\n";
                  mDynamicModelFile << "    g1=(double*)mxMalloc(" << block_triangular.ModelBlock->Block_List[i].Size*block_triangular.ModelBlock->Block_List[i].Size << "*sizeof(double));\n";
                  mDynamicModelFile << "    r=(double*)mxMalloc(" << block_triangular.ModelBlock->Block_List[i].Size << "*sizeof(double));\n";
                  if (!block_triangular.ModelBlock->Block_List[i].is_linear)
                    {
                      mDynamicModelFile << "    cvg=false;\n";
                      mDynamicModelFile << "    iter=0;\n";
                      mDynamicModelFile << "    while(!((cvg)||(iter>maxit_)))\n";
                      mDynamicModelFile << "      {\n";
                      mDynamicModelFile << "        res2=0;\n";
                      mDynamicModelFile << "        res1=0;\n";
                      mDynamicModelFile << "        max_res=0;\n";
                      mDynamicModelFile << "        for(it_=y_kmin;it_<periods+y_kmin;it_++)\n";
                      mDynamicModelFile << "          {\n";
                      mDynamicModelFile << "            Per_u_=(it_-y_kmin)*" << block_triangular.ModelBlock->Block_List[i].IM_lead_lag[block_triangular.ModelBlock->Block_List[i].Max_Lag + block_triangular.ModelBlock->Block_List[i].Max_Lead].u_finish + 1 << ";\n";
                      mDynamicModelFile << "            Per_y_=it_*y_size;\n";
                      mDynamicModelFile << "            Dynamic" << i + 1 << "(y, x, r, g1, g2);\n";
                      mDynamicModelFile << "            for(i=0;i<" << block_triangular.ModelBlock->Block_List[i].Size << ";i++)\n";
                      mDynamicModelFile << "              {\n";
                      mDynamicModelFile << "                if (max_res<fabs(r[i]))\n";
                      mDynamicModelFile << "                  max_res=fabs(r[i]);\n";
                      mDynamicModelFile << "                res2+=r[i]*r[i];\n";
                      mDynamicModelFile << "                res1+=fabs(r[i]);\n";
                      mDynamicModelFile << "              }\n";
                      mDynamicModelFile << "          }\n";
                      mDynamicModelFile << "        cvg=(max_res<solve_tolf);\n";
                      mDynamicModelFile << "        simulate(" << i << ", " << /*mod_param.endo_nbr*/symbol_table.endo_nbr << ", periods, y_kmin, y_kmax);\n";
                      mDynamicModelFile << "        iter++;\n";
                      mDynamicModelFile << "      }\n";
                      mDynamicModelFile << "    if (!cvg)\n";
                      mDynamicModelFile << "      {\n";
                      mDynamicModelFile << "        mexPrintf(\"Convergence not achieved in block " << i << ", after %d iterations\\n\",iter);\n";
                      mDynamicModelFile << "        mexErrMsgTxt(\"End of simulate\");\n";
                      mDynamicModelFile << "      }\n";
                    }
                  else
                    {
                      mDynamicModelFile << "    for(it_=y_kmin;it_<periods+y_kmin;it_++)\n";
                      mDynamicModelFile << "      {\n";
                      mDynamicModelFile << "        Per_u_=(it_-y_kmin)*" << block_triangular.ModelBlock->Block_List[i].IM_lead_lag[block_triangular.ModelBlock->Block_List[i].Max_Lag + block_triangular.ModelBlock->Block_List[i].Max_Lead].u_finish + 1 << ";\n";
                      mDynamicModelFile << "        Per_y_=it_*y_size;\n";
                      mDynamicModelFile << "        Dynamic" << i + 1 << "(y, x, r, g1, g2);\n";
                      mDynamicModelFile << "#ifdef PRINT_OUT\n";
                      mDynamicModelFile << "        for(j=0;j<" << block_triangular.ModelBlock->Block_List[i].IM_lead_lag[block_triangular.ModelBlock->Block_List[i].Max_Lag + block_triangular.ModelBlock->Block_List[i].Max_Lead].u_finish + 1 << ";j++)\n";
                      mDynamicModelFile << "          {\n";
                      mDynamicModelFile << "            mexPrintf(\" %f\",u[Per_u_+j]);\n";
                      mDynamicModelFile << "          }\n";
                      mDynamicModelFile << "        mexPrintf(\"\\n\");\n";
                      mDynamicModelFile << "#endif\n";
                      mDynamicModelFile << "      }\n";
                      mDynamicModelFile << "    simulate(" << i << ", " << /*mod_param.endo_nbr*/symbol_table.endo_nbr << ", it_, y_kmin, y_kmax);\n";
                    }
                  mDynamicModelFile << "#ifdef DEBUG\n";
                  mDynamicModelFile << "    for(it_=y_kmin;it_<periods+y_kmin;it_++)\n";
                  mDynamicModelFile << "      {\n";
                  mDynamicModelFile << "        for(i=0;i<Model_Block->List[" << i << "].Size;i++)\n";
                  mDynamicModelFile << "          {\n";
                  mDynamicModelFile << "            Per_y_=it_*y_size;\n";
                  mDynamicModelFile << "            mexPrintf(\" y[%d, %d]=%f \",it_,Model_Block->List[" << i << "].Variable[i],y[it_*" << /*mod_param.endo_nbr*/symbol_table.endo_nbr << "+Model_Block->List[" << i << "].Variable[i]]);\n";
                  mDynamicModelFile << "          }\n";
                  mDynamicModelFile << "        mexPrintf(\" \\n \");\n";
                  mDynamicModelFile << "      }\n";
                  mDynamicModelFile << "#endif\n";
                  mDynamicModelFile << "    mxFree(g1);\n";
                  mDynamicModelFile << "    mxFree(r);\n";
                  mDynamicModelFile << "    mxFree(u);\n";
                  mDynamicModelFile << "    //mexErrMsgTxt(\"Exit from Dynare\");\n";
                }
              prev_Simulation_Type=k;
            }
          if (open_par)
            {
              mDynamicModelFile << "#endif\n";
              mDynamicModelFile << "      }\n";
            }
          mDynamicModelFile << " }\n";
        }
      if (offset == 2)
        {
          // Writing the gateway routine
          mDynamicModelFile << " int max(int a, int b)\n";
          mDynamicModelFile << " {\n";
          mDynamicModelFile << "   if (a>b) return(a); else return(b);\n";
          mDynamicModelFile << " }\n\n\n";
          mDynamicModelFile << "/* The gateway routine */\n";
          mDynamicModelFile << "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])\n";
          mDynamicModelFile << "{\n";
          mDynamicModelFile << "  tModel_Block *Model_Block;\n";
          mDynamicModelFile << "  mxArray *M_, *oo_, *options_;\n";
          mDynamicModelFile << "  int i, j, row_y, col_y, row_x, col_x, x_FieldNumber;\n";
          mDynamicModelFile << "  mxArray *x_FieldByNumber;\n";
          mDynamicModelFile << "  double * pind ;\n";
          mDynamicModelFile << "\n";
          mDynamicModelFile << "  /* Gets model parameters from global workspace of Matlab */\n";
          mDynamicModelFile << "  M_ = mexGetVariable(\"global\",\"M_\");\n";
          mDynamicModelFile << "  if (M_ == NULL )\n";
          mDynamicModelFile << "    {\n";
          mDynamicModelFile << "      mexPrintf(\"Global variable not found : \");\n";
          mDynamicModelFile << "      mexErrMsgTxt(\"M_ \\n\");\n";
          mDynamicModelFile << "    }\n";
          mDynamicModelFile << "  /* Gets variables and parameters from global workspace of Matlab */\n";
          mDynamicModelFile << "  oo_ = mexGetVariable(\"global\",\"oo_\");\n";
          mDynamicModelFile << "  if (oo_ == NULL )\n";
          mDynamicModelFile << "    {\n";
          mDynamicModelFile << "      mexPrintf(\"Global variable not found : \");\n";
          mDynamicModelFile << "      mexErrMsgTxt(\"oo_ \\n\");\n";
          mDynamicModelFile << "    }\n";
          mDynamicModelFile << "  options_ = mexGetVariable(\"global\",\"options_\");\n";
          mDynamicModelFile << "  if (options_ == NULL )\n";
          mDynamicModelFile << "    {\n";
          mDynamicModelFile << "      mexPrintf(\"Global variable not found : \");\n";
          mDynamicModelFile << "      mexErrMsgTxt(\"options_ \\n\");\n";
          mDynamicModelFile << "    }\n";
          mDynamicModelFile << "  params = mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,\"params\")));\n";
          mDynamicModelFile << "  y= mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,\"endo_simul\")));\n";
          mDynamicModelFile << "  row_y=mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,\"endo_simul\")));\n";
          mDynamicModelFile << "  x= mxGetPr(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,\"exo_simul\")));\n";
          mDynamicModelFile << "  row_x=mxGetM(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,\"exo_simul\")));\n";
          mDynamicModelFile << "  col_x=mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,\"exo_simul\")));\n";
          if (compiler==GCC_COMPILE)
            {
              mDynamicModelFile << "  y_kmin=int(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,\"maximum_lag\"))))));\n";
              mDynamicModelFile << "  y_kmax=int(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,\"maximum_lead\"))))));\n";
              mDynamicModelFile << "  y_decal=max(0,y_kmin-int(floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,\"maximum_endo_lag\")))))));\n";
              mDynamicModelFile << "  periods=int(floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,\"periods\"))))));\n";
              mDynamicModelFile << "  maxit_=int(floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,\"maxit_\"))))));\n";
              mDynamicModelFile << "  slowc=double(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,\"slowc\")))));\n";
            }
          else
            {
              mDynamicModelFile << "  y_kmin=(int)floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,\"maximum_lag\")))));\n";
              mDynamicModelFile << "  y_kmax=(int)floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_,\"maximum_lead\")))));\n";
              mDynamicModelFile << "  periods=(int)floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,\"periods\")))));\n";
              mDynamicModelFile << "  maxit_=(int)floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,\"maxit_\")))));\n";
              mDynamicModelFile << "  slowc=(int)floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,\"slowc\")))));\n";
            }
          mDynamicModelFile << "  col_y=mxGetN(mxGetFieldByNumber(oo_, 0, mxGetFieldNumber(oo_,\"endo_simul\")));;\n";
          mDynamicModelFile << "  if (col_y<row_x)\n";
          mDynamicModelFile << "    {\n";
          mDynamicModelFile << "      row_y=row_y/row_x;\n";
          mDynamicModelFile << "      col_y=row_x;\n";
          mDynamicModelFile << "    }\n";
          mDynamicModelFile << "  solve_tolf=*(mxGetPr(mxGetFieldByNumber(options_, 0, mxGetFieldNumber(options_,\"dynatol\"))));\n";
          mDynamicModelFile << "  y_size=row_y;\n";
          mDynamicModelFile << "  x_size=periods+y_kmin+y_kmax;\n";
          mDynamicModelFile << "#ifdef DEBUG\n";
          mDynamicModelFile << "  for(j=0;j<periods+y_kmin+y_kmax;j++)\n";
          mDynamicModelFile << "    {\n";
          mDynamicModelFile << "      for(i=0;i<row_y;i++)\n";
          mDynamicModelFile << "        mexPrintf(\"y[%d,%d]=%f \",j,i,y[j*y_size+i]);\n";
          mDynamicModelFile << "      mexPrintf(\"\\n\");\n";
          mDynamicModelFile << "    }\n";
          mDynamicModelFile << "    mexPrintf(\"\\n\");\n";
          mDynamicModelFile << "    mexPrintf(\"x=%x\\n\",x);\n";

          mDynamicModelFile << "  for(j=0;j<periods+y_kmin+y_kmax;j++)\n";
          mDynamicModelFile << "    {\n";
          mDynamicModelFile << "      for(i=0;i<col_x;i++)\n";
          mDynamicModelFile << "        mexPrintf(\"x[%d,%d]=%f \",j,i,x[i*x_size+j]);\n";
          mDynamicModelFile << "      mexPrintf(\"\\n\");\n";
          mDynamicModelFile << "    }\n";
          mDynamicModelFile << "    mexPrintf(\"x[1]=%f\\n\",x[1]);\n";
          mDynamicModelFile << "#endif\n";

          mDynamicModelFile << "  /* Gets it_ from global workspace of Matlab */\n";
          mDynamicModelFile << "  it_ = (int) floor(mxGetScalar(mexGetVariable(\"global\", \"it_\")))-1;\n";
          mDynamicModelFile << "  /* Call the C subroutines. */\n";
          mDynamicModelFile << "  t0= pctimer();\n";
          mDynamicModelFile << "  Dynamic_Init(Model_Block);\n";
          mDynamicModelFile << "  t1= pctimer();\n";
          mDynamicModelFile << "  mexPrintf(\"Simulation Time=%f milliseconds\\n\",1000*(t1-t0));\n";
          if (compiler==GCC_COMPILE  )
            {
              mDynamicModelFile << "  if (SaveCode.is_open())\n";
              mDynamicModelFile << "    SaveCode.close();\n";
            }
          else
            {
              mDynamicModelFile << "  if (SaveCode)\n";
              mDynamicModelFile << "    fclose(SaveCode);\n";
            }
          mDynamicModelFile << "  if (nlhs>0)\n";
          mDynamicModelFile << "    {\n";
          mDynamicModelFile << "      plhs[0] = mxCreateDoubleMatrix(row_y, col_y, mxREAL);\n";
          mDynamicModelFile << "      pind = mxGetPr(plhs[0]);\n";
          mDynamicModelFile << "      for(i=0;i<row_y*col_y;i++)\n";
          mDynamicModelFile << "        pind[i]=y[i];\n";
          mDynamicModelFile << "    }\n";
          mDynamicModelFile << "}\n";
        }
      mDynamicModelFile.close();
      if (printed)
        cout << "done\n";
    }
}

void
ModelTree::writeDynamicModel(ostream &DynamicOutput, Model_Block *ModelBlock) const
{
  ostringstream lsymetric;       // Used when writing symetric elements in Hessian
  ostringstream model_output;    // Used for storing model equations
  ostringstream jacobian_output; // Used for storing jacobian equations
  ostringstream hessian_output;  // Used for storing Hessian equations
  ostringstream third_derivatives_output;

  if (offset != 2)
    writeTemporaryTerms(model_output, true);

  writeLocalParameters(model_output, true);

  if (offset != 2)
    writeModelEquations(model_output, true);
  else
    writeModelEquationsOrdered(model_output, true,ModelBlock);

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
            d1->writeOutput(jacobian_output, true, temporary_terms, offset);
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
        d2->writeOutput(hessian_output, true, temporary_terms, offset);
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
        d3->writeOutput(third_derivatives_output, true, temporary_terms, offset);
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
      if (offset!=2)
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
      else
        DynamicOutput << model_output.str();
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

inline void
ModelTree::Evaluate_Jacobian()
{
  int i=0;
  bool *IM;
  int a_variable_lag=-9999;
  for(first_derivatives_type::iterator it = first_derivatives.begin();
      it != first_derivatives.end(); it++)
    {
      if (variable_table.getType(it->first.second) == eEndogenous)
        {
          NodeID Id = it->second;
          Id->Evaluate();
          interprete_.u1 = interprete_.Stack.top();
          interprete_.Stack.pop();
          int eq=it->first.first;
          int var=variable_table.getSymbolID(it->first.second);
          int k1=variable_table.getLag(it->first.second);
          if (a_variable_lag!=k1)
            {
              IM=block_triangular.bGet_IM(k1);
              a_variable_lag=k1;
            }
          if (IM[eq*symbol_table.endo_nbr+var] && (fabs(interprete_.u1)<interprete_.cutoff))
            {
              //cout << "the coefficient related to variable " << var << " with lag " << k1 << " in equation " << eq << " is equal to " << interprete_.u1 << " and is set to 0 in the incidence matrix\n";
              block_triangular.unfill_IM(eq, var, k1);
              i++;
            }
        }
    }
  if (i>0)
    cout << i << " elements in the incidence matrices are below the cutoff (" << interprete_.cutoff << ") and are discarded\n";
}

inline void
ModelTree::BlockLinear(Model_Block *ModelBlock)
{
  int i,j,l,m;
  for(j = 0;j < ModelBlock->Size;j++)
    {
      if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_BACKWARD_COMPLETE ||
          ModelBlock->Block_List[j].Simulation_Type==SOLVE_FOREWARD_COMPLETE)
        {
          for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[0].size;i++)
            {
              int eq=ModelBlock->Block_List[j].IM_lead_lag[0].Equ_Index[i];
              int var=ModelBlock->Block_List[j].IM_lead_lag[0].Var_Index[i];
              first_derivatives_type::const_iterator it=first_derivatives.find(make_pair(eq,variable_table.getmVariableSelector(var,0)));
              NodeID Id = it->second;
              Id->collectEndogenous(Id);
              if (Id->present_endogenous_size()>0)
                {
                  for(l=0;l<ModelBlock->Block_List[j].Size;l++)
                    {
                      if (Id->present_endogenous_find(ModelBlock->Block_List[j].Variable[l],0))
                        {
                          ModelBlock->Block_List[j].is_linear=false;
                          goto follow;
                        }
                    }
                }
            }
        }
      else if (ModelBlock->Block_List[j].Simulation_Type==SOLVE_TWO_BOUNDARIES_COMPLETE)
        {
          for(m=0;m<=ModelBlock->Block_List[j].Max_Lead+ModelBlock->Block_List[j].Max_Lag;m++)
            {
              int k1=m-ModelBlock->Block_List[j].Max_Lag;
              for(i=0;i<ModelBlock->Block_List[j].IM_lead_lag[m].size;i++)
                {
                  int eq=ModelBlock->Block_List[j].IM_lead_lag[m].Equ_Index[i];
                  int var=ModelBlock->Block_List[j].IM_lead_lag[m].Var_Index[i];
                  first_derivatives_type::const_iterator it=first_derivatives.find(make_pair(eq,variable_table.getmVariableSelector(var,k1)));
                  NodeID Id = it->second;
                  Id->collectEndogenous(Id);
                  if (Id->present_endogenous_size()>0)
                    {
                      for(l=0;l<ModelBlock->Block_List[j].Size;l++)
                        {
                          if (Id->present_endogenous_find(ModelBlock->Block_List[j].Variable[l],k1))
                            {
                              ModelBlock->Block_List[j].is_linear=false;
                              goto follow;
                            }
                        }
                    }
                }
            }
        }
    follow:
      i=0;
    }
}

void
ModelTree::computingPass()
{
  cout << equations.size() << " equation(s) found" << endl;

  // Sorting variable table
  variable_table.Sort();

  variable_table.setmVariableSelector();

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

  if (offset == 2)
    {
      int Size;
      int HSize;
      int *Table=variable_table.GetVariableTable(&Size,&HSize);

      interprete_.create_id_map(Table,Size,HSize);
      Evaluate_Jacobian();

      if (block_triangular.bt_verbose)
        {
          cout << "The gross incidence matrix \n";
          block_triangular.Print_IM( symbol_table.endo_nbr);
        }
      block_triangular.SetVariableTable(Table, Size, HSize);
      block_triangular.Normalize_and_BlockDecompose_Static_0_Model();
      BlockLinear(block_triangular.ModelBlock);

      computeTemporaryTermsOrdered(order, block_triangular.ModelBlock);
    }
  else
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
  if (offset == 1)
    writeDynamicMFile(basename + "_dynamic");
  else
    writeDynamicCFile(basename + "_dynamic");
}
