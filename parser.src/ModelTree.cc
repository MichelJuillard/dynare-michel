/*! \file
  \version 1.0
  \date 04/09/2004
  \par This file implements the ModelTree class methodes.
*/
//------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <numeric>
#include <stdio.h>
#include <map>
#include <time.h>
using namespace std;
//------------------------------------------------------------------------------
#include "VariableTable.hh"
#include "DynareBison.hh"
#include "NumericalConstants.hh"
#include "ModelTree.hh"
#include "ModelParameters.hh"
#include "Interface.hh"
//------------------------------------------------------------------------------
ostringstream ModelTree::output;
//------------------------------------------------------------------------------

inline NodeID MetaToken::getDerivativeAddress(int iVarID)
{
  std::map<int, NodeID, std::less<int> >::iterator iter = d1.find(iVarID);
  if (iter == d1.end())
    // No entry in map, derivative is therefore null
    if (op_code == EQUAL)
      return DataTree::ZeroEqZero;
    else
      return DataTree::Zero;
  else
    return iter->second;
}

ModelTree::ModelTree()
{
  computeJacobian = false;
  computeJacobianExo = false;
  computeHessian = false;
}

//------------------------------------------------------------------------------
ModelTree::~ModelTree()
{
  // Empty
}

//------------------------------------------------------------------------------
void ModelTree::OpenMFiles(string iModelFileName1, string iModelFileName2)
{
  if (iModelFileName1.size())
    {
      iModelFileName1 += interfaces::function_file_extension();
      mStaticModelFile.open(iModelFileName1.c_str(),ios::out|ios::binary);
      if (!mStaticModelFile.is_open())
        {
          cout << "ModelTree::Open : Error : Can't open file " << iModelFileName1
               << " for writing\n";
          exit(-1);
        }
      iModelFileName1.erase(iModelFileName1.end()-2,iModelFileName1.end());
      //Writing comments and function definition command
      mStaticModelFile << "function [residual, g1] = " <<  iModelFileName1 << "( y, x )\n";
      mStaticModelFile << interfaces::comment()+"\n"+interfaces::comment();
      mStaticModelFile << "Status : Computes static model for Dynare\n" << interfaces::comment() << "\n";
      mStaticModelFile << interfaces::comment();
      mStaticModelFile << "Warning : this file is generated automatically by Dynare\n";
      mStaticModelFile << interfaces::comment();
      mStaticModelFile << "  from model file (.mod)\n\n";
      if (iModelFileName2.size() && (computeJacobian||computeJacobianExo||computeHessian))
        {
          iModelFileName2 += interfaces::function_file_extension();
          mDynamicModelFile.open(iModelFileName2.c_str(),ios::out|ios::binary);
          if (!mDynamicModelFile.is_open())
            {
              cout << "ModelTree::Open : Error : Can't open file " << iModelFileName2
                   << " for writing\n";
              exit(-1);
            }
          iModelFileName2.erase(iModelFileName2.end()-2,iModelFileName2.end());
          mDynamicModelFile << "function [residual, g1, g2] = " <<  iModelFileName2 << "(y, x)\n";
          mDynamicModelFile << interfaces::comment()+"\n"+interfaces::comment();
          mDynamicModelFile << "Status : Computes dynamic model for Dynare\n" << interfaces::comment() << "\n";
          mDynamicModelFile << interfaces::comment();
          mDynamicModelFile << "Warning : this file is generated automatically by Dynare\n";
          mDynamicModelFile << interfaces::comment();
          mDynamicModelFile << "  from model file (.mod)\n\n";

        }
    }
  else
    {
      cout << "ModelTree::Open : Error : Missing file name\n";
      exit(-1);
    }
}

//------------------------------------------------------------------------------
void ModelTree::OpenCFiles(string iModelFileName1, string iModelFileName2)
{
  if (iModelFileName1.size())
    {
      iModelFileName1 += ".c";
      mStaticModelFile.open(iModelFileName1.c_str(),ios::out|ios::binary);
      if (!mStaticModelFile.is_open())
        {
          cout << "ModelTree::Open : Error : Can't open file " << iModelFileName1
               << " for writing\n";
          exit(-1);
        }
      iModelFileName1.erase(iModelFileName1.end()-2,iModelFileName1.end());
      mStaticModelFile << "/*\n";
      mStaticModelFile << " *" << iModelFileName1 << ".c  : Computes static model for Dynare\n";
      mStaticModelFile << " * Warning : this file is generated automatically by Dynare\n";
      mStaticModelFile << " *           from model file (.mod)\n\n";
      mStaticModelFile << " */\n";
      mStaticModelFile << "#include <math.h>\n";
      mStaticModelFile << "#include \"mex.h\"\n";
      // A flobal variable for model parameters
      mStaticModelFile << "double *params;\n";
      if (iModelFileName2.size() && (computeJacobian||computeJacobianExo||computeHessian))
        {
          iModelFileName2 += ".c";
          mDynamicModelFile.open(iModelFileName2.c_str(),ios::out|ios::binary);
          if (!mDynamicModelFile.is_open())
            {
              cout << "ModelTree::Open : Error : Can't open file " << iModelFileName2
                   << " for writing\n";
              exit(-1);
            }
          iModelFileName2.erase(iModelFileName2.end()-2,iModelFileName2.end());
          mDynamicModelFile << "/*\n";
          mDynamicModelFile << " *" << iModelFileName2 << ".c  : Computes dynamic model for Dynare\n";
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
        }
    }
  else
    {
      cout << "ModelTree::Open : Error : Missing file name\n";
      exit(-1);
    }
}

//------------------------------------------------------------------------------
void ModelTree::SaveMFiles()
{
  if (mStaticModelFile.is_open())
    {
      mStaticModelFile << StaticOutput.str();
      interfaces::function_close();
      mStaticModelFile.close();
    }
  if (mDynamicModelFile.is_open() && (computeJacobian||computeJacobianExo||computeHessian))
    {
      mDynamicModelFile << DynamicOutput.str();
      interfaces::function_close();
      mDynamicModelFile.close();
    }
}

//------------------------------------------------------------------------------
void ModelTree::SaveCFiles()
{
  if (mStaticModelFile.is_open())
    {
      // Writing the function Static
      mStaticModelFile << StaticOutput.str();
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
      mStaticModelFile << "      plhs[0] = mxCreateDoubleMatrix(" << ModelParameters::eq_nbr << ",1, mxREAL);\n";
      mStaticModelFile << "     /* Create a C pointer to a copy of the output matrix residual. */\n";
      mStaticModelFile << "     residual = mxGetPr(plhs[0]);\n";
      mStaticModelFile << "  }\n\n";
      mStaticModelFile << "  g1 = NULL;\n";
      mStaticModelFile << "  if (nlhs >= 2)\n";
      mStaticModelFile << "  {\n";
      mStaticModelFile << "      /* Set the output pointer to the output matrix g1. */\n";
      mStaticModelFile << "      plhs[1] = mxCreateDoubleMatrix(" << ModelParameters::eq_nbr << ", " << ModelParameters::endo_nbr << ", mxREAL);\n";
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
  if (mDynamicModelFile.is_open() && (computeJacobian||computeJacobianExo||computeHessian))
    {
      // Writing the function body
      mDynamicModelFile << DynamicOutput.str();
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
      mDynamicModelFile << "     plhs[0] = mxCreateDoubleMatrix(" << ModelParameters::eq_nbr << ",1, mxREAL);\n";
      mDynamicModelFile << "     /* Create a C pointer to a copy of the output matrix residual. */\n";
      mDynamicModelFile << "     residual = mxGetPr(plhs[0]);\n";
      mDynamicModelFile << "  }\n\n";
      mDynamicModelFile << "  g1 = NULL;\n";
      mDynamicModelFile << "  if (nlhs >= 2)\n";
      mDynamicModelFile << "  {\n";
      mDynamicModelFile << "     /* Set the output pointer to the output matrix g1. */\n";
      if (computeJacobianExo)
        mDynamicModelFile << "     plhs[1] = mxCreateDoubleMatrix(" << ModelParameters::eq_nbr << ", " << VariableTable::size() << ", mxREAL);\n";
      else if (computeJacobian)
        mDynamicModelFile << "     plhs[1] = mxCreateDoubleMatrix(" << ModelParameters::eq_nbr << ", " << ModelParameters::var_endo_nbr << ", mxREAL);\n";
      mDynamicModelFile << "     /* Create a C pointer to a copy of the output matrix g1. */\n";
      mDynamicModelFile << "     g1 = mxGetPr(plhs[1]);\n";
      mDynamicModelFile << "  }\n\n";
      mDynamicModelFile << "  g2 = NULL;\n";
      mDynamicModelFile << " if (nlhs >= 3)\n";
      mDynamicModelFile << "  {\n";
      mDynamicModelFile << "     /* Set the output pointer to the output matrix g2. */\n";
      mDynamicModelFile << "     plhs[2] = mxCreateDoubleMatrix(" << ModelParameters::eq_nbr << ", " << VariableTable::size()*VariableTable::size() << ", mxREAL);\n";
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
}

//------------------------------------------------------------------------------

void ModelTree::derive(int iOrder)
{
  NodeID  lToken;                // To store current working token
  NodeID    lD1, lD2;            // To store derivative arguments of current argument
  NodeID   lArg1, lArg2;         // To store current arguments
  Type    lType1;                // Type of first argument
  NodeID   t1,t11,t12,t13,
    t14, t15;                    // To store temoporary result arguments
  TreeIterator  BeginIT;         // Iterator of the 1st token to derive
  TreeIterator  EndIT;           // Iterator of the last token to derive
  TreeIterator  currentIT;       // Iterator (counter) for model tree loop

  vector<NodeID> EqualTokenIDs;  // IDs of "equal token" in model Tree
  // Capturing equation IDs
  for (currentIT = BeginModel; currentIT != mModelTree.end(); currentIT++)
    {
      if ((*currentIT)->op_code == EQUAL)
        {
          EqualTokenIDs.push_back(*currentIT);
          // Equation is forced to be in Model Tree as referenced
          // This is useful to remove symmetric elements
          (*currentIT)->reference_count[0]++;
        }
    }
  //std::cout << "size " << EqualTokenIDs.size() << "\n";
  mDerivativeIndex.resize(iOrder);
  // Uncomment this to print model tree data
  /*
  //cout << "ModelTree==================================\n";
  for (currentIT = mModelTree.begin(); currentIT != mModelTree.end(); currentIT++)
  {
  lToken = *currentIT;
  int ID = lToken->idx;
  cout << ID << ":" << lToken << "->" << lToken->id1 << " " << lToken->type1 << " " <<
  lToken->id2 << " " << lToken->op_code << "\n";
  }
  */

  EndIT = mModelTree.begin();
  EndIT--;
  cout << "Processing derivation ...\n";
  // loop on order of derivation
  for(int Order = 1; Order <= iOrder; Order++)
    {

      cout << "\tProcessing Order " << Order << "... ";
      current_order = Order;
      BeginIT = EndIT;
      BeginIT++;
      EndIT = mModelTree.end();
      EndIT--;
      // Adding a reference counter for current order to tokens in mModelTree
      // and updating them

      for (TreeIterator it = mModelTree.begin(); it != mModelTree.end(); it++)
        {
          int s = (*it)->reference_count.size();
          for (int i = s; i <= current_order; i++)
            {
              int rc = (*it)->reference_count[i-1];
              (*it)->reference_count.push_back(rc);
            }
        }

      // Loop on tokens
      for (currentIT = BeginIT;; currentIT++)
        {
          lToken = *currentIT;

          // Loop on variables of derivation which may give non null result for this token
          for (vector<int>::iterator varIt = lToken->non_null_derivatives.begin();
               varIt != lToken->non_null_derivatives.end(); varIt++)
            {
              int var = *varIt;
              //cout << "Token " << (*currentIT)->idx << endl;
              if (accumulate((*currentIT)->reference_count.begin(), (*currentIT)->reference_count.end(),0) > 0)
                {
                  lArg1 = lToken->id1;
                  lArg2 = lToken->id2;
                  lType1 = lToken->type1;
                  lD1 = DeriveArgument(lArg1, lType1, var);
                  lD2 = Zero;
                  if (lArg2 != NullID)
                    lD2 = DeriveArgument(lArg2, eTempResult, var);
                  // Case where token is a terminal
                  if (lToken->op_code == NoOpCode)
                    {
                      (*currentIT)->setDerivativeAddress(lD1, var);
                    }
                  else if (lD1 == Zero && lD2 == Zero)
                    {
                      (*currentIT)->setDerivativeAddress(Zero, var);
                    }
                  else
                    {
                      switch (lToken->op_code)
                        {
                        case UMINUS:
                          t1 = AddUMinus(lD1);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case PLUS:
                          t1 = AddPlus(lD1, lD2);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case MINUS:
                          t1 = AddMinus(lD1, lD2);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case TIMES:
                          t11 = AddTimes(lD1, lArg2);
                          t12 = AddTimes(lD2, lArg1);
                          t1 = AddPlus(t11, t12);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case DIVIDE:
                          t11 = AddTimes(lD1, lArg2);
                          t12 = AddTimes(lD2, lArg1);
                          t13 = AddMinus(t11, t12);
                          t14 =  AddTimes(lArg2, lArg2);
                          t1 = AddDivide(t13, t14);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case SQRT:
                          t11 = AddPlus(*currentIT, *currentIT);
                          t1 = AddDivide(lD1, t11);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case POWER:
                          if (lD2 == Zero)
                            {
                              if (lD1 == Zero)
                                t1 = Zero;
                              else
                                {
                                  t11 = AddMinus(lArg2, One);
                                  t12 = AddPower(lArg1, t11);
                                  t13 = AddTimes(lArg2, t12);
                                  t1 = AddTimes(lD1, t13);
                                }
                            }
                          else
                            {
                              t11 = AddLog(lArg1);
                              t12 = AddTimes(lD2, t11);
                              t13 = AddTimes(lD1, lArg2);
                              t14 =  AddDivide(t13, lArg1);
                              t15 = AddPlus(t12, t14);
                              t1 = AddTimes(t15, *currentIT);
                            }
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case EXP:
                          t1 = AddTimes(lD1, *currentIT);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case LOG:
                          t1 = AddDivide(lD1, lArg1);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case LOG10:
                          t11 = AddExp(One);
                          t12 = AddLog10(t11);
                          t13 = AddDivide(lD1, lArg1);
                          t1 = AddTimes(t12, t13);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case COS:
                          t11 = AddSin(lArg1);
                          t12 = AddUMinus(t11);
                          t1 = AddTimes( lD1, t12);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case SIN:
                          t11 = AddCos(lArg1);
                          t1 = AddTimes(lD1,t11);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case TAN:
                          t11 = AddTimes(*currentIT, *currentIT);
                          t12 = AddPlus(t11, One);
                          t1 = AddTimes(lD1, t12);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case ACOS:
                          t11 = AddSin(*currentIT);
                          t12 = AddDivide(lD1, t11);
                          t1 = AddUMinus(t12);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case ASIN:
                          t11 = AddCos(*currentIT);
                          t1 = AddDivide(lD1, t11);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case ATAN:
                          t11 = AddTimes(lArg1, lArg1);
                          t12 = AddPlus(One, t11);
                          t1 = AddDivide(lD1, t12);
                          (*currentIT)->setDerivativeAddress(t1,var);
                          break;
                        case COSH:
                          t11 = AddSinH(lArg1);
                          t1 = AddTimes( lD1,t11);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case SINH:
                          t11 = AddCosH(lArg1);
                          t1 = AddTimes( lD1, t11);
                          (*currentIT)->setDerivativeAddress(t1, var);
                          break;
                        case TANH:
                          t11 = AddTimes(*currentIT, *currentIT);
                          t12 = AddMinus(One, t11);
                          t1 = AddTimes(lD1, t12);
                          (*currentIT)->setDerivativeAddress(t1,var);
                          break;
                        case ACOSH:
                          t11 = AddSinH(*currentIT);
                          t1 = AddDivide(lD1, t11);
                          (*currentIT)->setDerivativeAddress(t1,var);
                          break;
                        case ASINH:
                          t11 = AddCosH(*currentIT);
                          t1 = AddDivide(lD1, t11);
                          (*currentIT)->setDerivativeAddress(t1,var);
                          break;
                        case ATANH:
                          t11 = AddTimes(lArg1, lArg1);
                          t12 = AddMinus(One, t11);
                          t1 = AddTimes(lD1, t12);
                          (*currentIT)->setDerivativeAddress(t1,var);
                          break;
                        case EQUAL:
                          // Force the derivative to have zero on right hand side
                          // (required for setStaticModel and setDynamicModel)
                          if (lD1 == Zero && lD2 != Zero)
                            {
                              t11 = AddUMinus(lD2);
                              t1 = AddEqual(t11, Zero);
                            }
                          else if (lD1 != Zero && lD2 == Zero)
                            t1 = AddEqual(lD1, Zero);
                          else
                            {
                              t11 = AddMinus(lD1, lD2);
                              t1 = AddEqual(t11, Zero);
                            }
                          // The derivative is forced to be in Model Tree as referenced
                          // This is useful to remove symmetric elements
                          IncrementReferenceCount(t1);
                          lToken->setDerivativeAddress(t1, var);
                          break;
                        }
                    }
                }
            }
          if (currentIT == EndIT)
            break;
        }

      // Filling mDerivativeIndex
      // Loop on variables of derivation, skipping symmetric elements
      for (int var = 0; var < VariableTable::size(); var++)
        {
          int starti = var*Order*(Order-1)*ModelParameters::eq_nbr/2;
          for (unsigned int i = starti; i < EqualTokenIDs.size() ; i++ )
            {
              t1 = EqualTokenIDs[i]->getDerivativeAddress(var);
              if (Order == 1)
                mDerivativeIndex[0].push_back(DerivativeIndex(t1, i-starti, var));
              else if (Order == 2)
                {
                  int var1 = VariableTable::getSortID(i/ModelParameters::eq_nbr);
                  int var2 = VariableTable::getSortID(var);
                  mDerivativeIndex[1].push_back(DerivativeIndex(t1,
                                                                i-ModelParameters::eq_nbr*(i/ModelParameters::eq_nbr),
                                                                var1*VariableTable::size()+var2));
                }
            }
        }

      // Uncomment to debug : prints unreferenced tokens
      /*
        cout << "Order : " << Order << "\n";
        for (TokenCount = BeginModel; TokenCount < mModelTree.size() ; TokenCount++ )
        {
        if (accumulate(mModelTree[TokenCount].reference_count.begin(),mModelTree[TokenCount].reference_count.end(),0) == 0)
        cout << "\tNot referenced : token ID :" << TokenCount << endl;
        }
      */
      // Uncomment this to debug : mDerivative(1and2)Index data
      // before removing unreferenced tokens
      /*
        cout << "Contenence of mDerivative1Index\n";
        for (int i=0; i< mDerivativeIndex[0].size();i++)
        //if (mDerivativeIndex[0][i].token_id != 3)
        cout << "\t" << mDerivativeIndex[0][i].token_id << endl;
        cout << "Contenence of mDerivative2Index\n";
        for (int i=0; i< mDerivativeIndex[1].size();i++)
        //if (mDerivativeIndex[1][i].token_id != 3)
        cout << "\t" << mDerivativeIndex[1][i].token_id << endl;
      */
      //cout << "Removing unreferenced tokens range ids :" << CurrentID << " - " << mModelTree.size()-1 << endl;
      // Removing unreferenced tokens in last derivative
      // RemoveUnref(CurrentID, mModelTree.size()-1, Order);
      // Decrementing reference couter of unreferenced tokens in last derivative
      //DecrementUnref(CurrentID, mModelTree.size()-1, Order);
      /*
        cout << "Order : " << Order << "\n";
        for (TokenCount = BeginModel; TokenCount < mModelTree.size() ; TokenCount++ )
        {
        if (accumulate(mModelTree[TokenCount].reference_count.begin(),mModelTree[TokenCount].reference_count.end(),0) == 0)
        cout << "\tNot referenced : token ID :" << TokenCount << endl;
        }
      */
      EqualTokenIDs.clear();
      // Updating EqualTokenIDs
      for (unsigned int i=0; i< mDerivativeIndex[Order-1].size();i++)
        {
          EqualTokenIDs.push_back(mDerivativeIndex[Order-1][i].token_id);
        }

      // Uncomment this to debug : mDerivative(1and2)Index data
      // after removing unreferenced tokens
      /*
        cout << "Contenence of mDerivative1Index after removing\n";
        for (int i=0; i< mDerivativeIndex[0].size();i++)
        if (mDerivativeIndex[0][i].token_id != 3)
        cout << "\t" << mDerivativeIndex[0][i].token_id << endl;
        cout << "Contenence of mDerivative2Index after removing\n";
        for (int i=0; i< mDerivativeIndex[1].size();i++)
        //if (mDerivativeIndex[1][i].token_id != 3)
        cout << "\t" << mDerivativeIndex[1][i].token_id << endl;
      */
      cout << "done" << endl;
    }
}

//------------------------------------------------------------------------------
inline NodeID ModelTree::DeriveArgument(NodeID iArg, Type iType, int iVarID)
{
  switch(iType)
    {
    case eTempResult      :
      return iArg->getDerivativeAddress(iVarID);
    case eExogenous       :
    case eExogenousDet      :
    case eEndogenous      :
    case eRecursiveVariable   :
      if ((long int) iArg == iVarID)
        //if ((VariableTable::getSymbolID(iArg) == VariableTable::getSymbolID(iVarID)) &&
        //   (VariableTable::getType(iArg) == VariableTable::getType(iVarID)))
        {
          /*
            cout << SymbolTable::getNameByID(iType,
            VariableTable::getSymbolID(iArg)) << endl;
            cout << SymbolTable::getNameByID(iType,
            VariableTable::getSymbolID(iVarID)) << endl;
          */
          return One;
        }
      else
        {
          return Zero;
        }
      break;
    case eNumericalConstant   :
    case eParameter       :
      return Zero;
    case eLocalParameter      :
      return Zero;
    case eUNDEF         :
      return NullID;
    default       :
      cout << "ModelTree::DeriveArgument : Error: Unknown Type\n";
      exit(-1);
    };

}

//------------------------------------------------------------------------------
string  ModelTree::setStaticModel(void)
{
  TreeIterator tree_it;
  int lEquationNBR = 0;
  ostringstream model_output;    // Used for storing model equations
  ostringstream model_tmp_output;// Used for storing tmp expressions for model equations
  ostringstream jacobian_output; // Used for storing jacobian equations
                                 // Used for storing tmp expressions for jacobian equations
  ostringstream jacobian_tmp_output;

  int d = current_order;         // Minimum number of times a temparary expression apears in equations
  int EquationNBR;               // Number of model equations

  EquationNBR = ModelParameters::eq_nbr;
  // Reference count of token "0=0" is set to 0
  // Not to be printed as a temp expression
  fill(ZeroEqZero->reference_count.begin(),
       ZeroEqZero->reference_count.end(),0);
  // Writing model Equations
  current_order = 1;
  tree_it = BeginModel;
  for (; tree_it != mModelTree.end(); tree_it++)
    {
      if ((*tree_it)->op_code == EQUAL || (*tree_it)->op_code == ASSIGN )
        {
          if ((*tree_it)->id1->type1 == eLocalParameter)
            {
              model_output << getExpression((*tree_it)->id1, eStaticEquations, lEquationNBR);
              model_output << " = ";
              model_output << getExpression((*tree_it)->id2, eStaticEquations, lEquationNBR) << ";" << endl;
            }
          else if  (lEquationNBR < ModelParameters::eq_nbr)
            {
              model_output << "lhs =";
              model_output << getExpression((*tree_it)->id1, eStaticEquations, lEquationNBR) << ";" << endl;
              model_output << "rhs =";
              model_output << getExpression((*tree_it)->id2, eStaticEquations, lEquationNBR) << ";" << endl;
              model_output << "residual" << lpar << lEquationNBR+1 << rpar << "= lhs-rhs;" << endl;
              lEquationNBR++;
            }
          else
            break;
        }
      else if ((*tree_it)->op_code != NoOpCode)
        {
          if (optimize(*tree_it))
            {
              model_output << "T" << (*tree_it)->idx << "=" << getExpression(*tree_it, eStaticEquations, lEquationNBR) << ";" << endl;
              (*tree_it)->tmp_status = 1;
            }
          else
            {
              (*tree_it)->tmp_status = 0;
            }
        }
    }

  // Writing Jacobian for endogenous variables without lag
  for(; tree_it != mModelTree.end(); tree_it++)
    {
      if ((*tree_it)->op_code != NoOpCode
          && (*tree_it)->op_code != EQUAL
          && (*tree_it)->op_code != ASSIGN)
        {
          if (optimize(*tree_it) == 1)
            {
              jacobian_output << "T" << (*tree_it)->idx << "=" << getExpression(*tree_it, eStaticEquations, lEquationNBR) << ";" << endl;
              (*tree_it)->tmp_status = 1;
            }
          else
            {
              (*tree_it)->tmp_status = 0;
            }
        }
    }

  lEquationNBR = 0;
  for (unsigned int i = 0; i < mDerivativeIndex[0].size(); i++)
    {
      if (VariableTable::getType(mDerivativeIndex[0][i].derivators) == eEndogenous)
        {
          NodeID startJacobian = mDerivativeIndex[0][i].token_id;
          if (startJacobian != ZeroEqZero)
            {
              string exp = getExpression(startJacobian->id1, eStaticDerivatives);
              ostringstream g1;
              g1 << "  g1" << lpar << mDerivativeIndex[0][i].equation_id+1 << ", " <<
                VariableTable::getSymbolID(mDerivativeIndex[0][i].derivators)+1 << rpar;
              jacobian_output << g1.str() << "=" <<  g1.str() << "+" << exp << ";\n";
            }
        }
    }

  // Writing ouputs
  if (offset == 1)
    {
      StaticOutput << "global M_ \n";
      StaticOutput << "if M_.param_nbr > 0\n  params = M_.params;\nend\n";

      StaticOutput << "  residual = zeros( " << ModelParameters::eq_nbr << ", 1);\n";
      StaticOutput << "\n\t"+interfaces::comment()+"\n\t"+interfaces::comment();
      StaticOutput << "Model equations\n\t";
      StaticOutput << interfaces::comment() + "\n\n";
      StaticOutput << model_output.str();
      StaticOutput << "if ~isreal(residual)\n";
      StaticOutput << "  residual = real(residual)+imag(residual).^2;\n";
      StaticOutput << "end\n";
      StaticOutput << "if nargout >= 2,\n";
      StaticOutput << "  g1 = " <<
        "zeros(" << ModelParameters::eq_nbr << ", " <<
        ModelParameters::endo_nbr << ");\n" ;
      StaticOutput << "\n\t"+interfaces::comment()+"\n\t"+interfaces::comment();
      StaticOutput << "Jacobian matrix\n\t";
      StaticOutput << interfaces::comment() + "\n\n";
      StaticOutput << jacobian_output.str();
      StaticOutput << "  if ~isreal(g1)\n";
      StaticOutput << "    g1 = real(g1)+2*imag(g1);\n";
      StaticOutput << "  end\n";
      StaticOutput << "end\n";
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
  current_order = d;
  return StaticOutput.str();
}

//------------------------------------------------------------------------------
string  ModelTree::setDynamicModel(void)
{
  TreeIterator tree_it;
  int lEquationNBR = 0;
  ostringstream lsymetric;       // Used when writing symetric elements
  ostringstream model_output;    // Used for storing model equations
  ostringstream model_tmp_output;// Used for storing tmp expressions for model equations
  ostringstream jacobian_output; // Used for storing jacobian equations
                                 // Used for storing tmp expressions for jacobian equations
  ostringstream jacobian_tmp_output;
  ostringstream hessian_output;  // Used for storing Hessian equations
                                 // Used for storing tmp expressions for Hessian equations
  ostringstream hessian_tmp_output;

  int d = current_order;

  // Reference count of token "0=0" is set to 0
  // Not to be printed as a temp expression
  fill(ZeroEqZero->reference_count.begin(),
       ZeroEqZero->reference_count.end(),0);

  // Clearing output string
  model_output.str("");
  jacobian_output.str("");
  // Getting equations from model tree
  // Starting from the end of equation
  // Searching for the next '=' operator,
  current_order = 1;
  lEquationNBR = 0;
  cout << "\tequations .. ";
  tree_it = BeginModel;
  for (; tree_it != mModelTree.end(); tree_it++)
    {
      if ((*tree_it)->op_code == EQUAL || (*tree_it)->op_code == ASSIGN)
        {
          if ((*tree_it)->id1->type1 == eLocalParameter)
            {
              model_output << getExpression((*tree_it)->id1, eStaticEquations, lEquationNBR);
              model_output << " = ";
              model_output << getExpression((*tree_it)->id2, eStaticEquations, lEquationNBR) << ";" << endl;
            }
          else if  (lEquationNBR < ModelParameters::eq_nbr)
            {
              model_output << "lhs =";
              model_output << getExpression(((*tree_it)->id1), eDynamicEquations, lEquationNBR) << ";" << endl;
              model_output << "rhs =";
              model_output << getExpression(((*tree_it)->id2), eDynamicEquations, lEquationNBR) << ";" << endl;
              model_output << "residual" << lpar << lEquationNBR+1 << rpar << "= lhs-rhs;" << endl;
              lEquationNBR++;
            }
          else
            break;
        }
      else if ((*tree_it)->op_code != NoOpCode)
        {
          if (optimize(*tree_it))
            {
              model_output << "T" << (*tree_it)->idx << "=" << getExpression(*tree_it, eDynamicEquations, lEquationNBR) << ";" << endl;
              (*tree_it)->tmp_status = 1;
            }
          else
            {
              (*tree_it)->tmp_status = 0;
            }
        }
    }

  for(; tree_it != mModelTree.end(); tree_it++)
    {
      if ((*tree_it)->op_code != NoOpCode
          && (*tree_it)->op_code != EQUAL
          && (*tree_it)->op_code != ASSIGN)
        {
          if (optimize(*tree_it) == 1)
            {
              jacobian_output << "T" << (*tree_it)->idx << "=" << getExpression(*tree_it, eDynamicEquations, lEquationNBR) << ";" << endl;
              (*tree_it)->tmp_status = 1;
            }
          else
            {
              (*tree_it)->tmp_status = 0;
            }
        }
    }
  cout << "done \n";
  // Getting Jacobian from model tree
  if (computeJacobian || computeJacobianExo)
    {
      cout << "\tJacobian .. ";
      for(; tree_it != mModelTree.end(); tree_it++)
        {
          if ((*tree_it)->op_code != NoOpCode
              && (*tree_it)->op_code != EQUAL
              && (*tree_it)->op_code != ASSIGN)
            {
              if (optimize(*tree_it) == 1)
                {
                  jacobian_output << "T" << (*tree_it)->idx << "=" << getExpression(*tree_it, eDynamicEquations, lEquationNBR) << ";" << endl;
                  (*tree_it)->tmp_status = 1;
                }
              else
                {
                  (*tree_it)->tmp_status = 0;
                }
            }
        }

      lEquationNBR = 0;
      for (unsigned int i = 0; i < mDerivativeIndex[0].size(); i++)
        {
          if (computeJacobianExo || VariableTable::getType(mDerivativeIndex[0][i].derivators) == eEndogenous)
            {
              NodeID startJacobian = mDerivativeIndex[0][i].token_id;
              if (startJacobian != ZeroEqZero)
                {
                  string exp = getExpression(startJacobian->id1, eDynamicDerivatives);
                  ostringstream g1;
                  g1 << "  g1" << lpar << mDerivativeIndex[0][i].equation_id+1 << ", " <<
                    VariableTable::getSortID(mDerivativeIndex[0][i].derivators)+1 << rpar;
                  jacobian_output << g1.str() << "=" <<  g1.str() << "+" << exp << ";\n";
                }
            }
        }
      cout << "done \n";
    }

  if (computeHessian)
    {
      // Getting Hessian from model tree
      // Starting from the end of equation
      // Searching for the next '=' operator,
      lEquationNBR = 0;
      cout << "\tHessian .. ";
      for (unsigned int i = 0; i < mDerivativeIndex[1].size(); i++)
        {
          NodeID startHessian = mDerivativeIndex[1][i].token_id;
          //cout << "ID = " << startHessian << " exp = " << exp << "\n";
          if (startHessian != ZeroEqZero)
            {
              string exp = getExpression(startHessian->id1, eDynamicDerivatives);

              int varID1 = mDerivativeIndex[1][i].derivators/VariableTable::size();
              int varID2 = mDerivativeIndex[1][i].derivators-varID1*VariableTable::size();
              hessian_output << "  g2" << lpar << mDerivativeIndex[1][i].equation_id+1 << ", " <<
                mDerivativeIndex[1][i].derivators+1 << rpar << " = " << exp << ";\n";
              // Treating symetric elements
              if (varID1 != varID2)
                lsymetric <<  "  g2" << lpar << mDerivativeIndex[1][i].equation_id+1 << ", " <<
                  varID2*VariableTable::size()+varID1+1 << rpar << " = " <<
                  "g2" << lpar << mDerivativeIndex[1][i].equation_id+1 << ", " <<
                  mDerivativeIndex[1][i].derivators+1 << rpar << ";\n";
            }

        }
      cout << "done \n";
    }
  int nrows = ModelParameters::eq_nbr;
  int nvars = ModelParameters::var_endo_nbr+ModelParameters::exo_nbr+ModelParameters::exo_det_nbr;
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
          int ncols = nvars*nvars;
          DynamicOutput << "  g2 = " <<
            "sparse([],[],[]," << nrows << ", " << ncols << ", " <<
            5*ncols << ");\n";
          DynamicOutput << "\n\t"+interfaces::comment()+"\n\t"+interfaces::comment();
          DynamicOutput << "Hessian matrix\n\t";
          DynamicOutput << interfaces::comment() + "\n\n";
          DynamicOutput << hessian_output.str() << lsymetric.str();
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
  current_order = d;
  return DynamicOutput.str();
}

//------------------------------------------------------------------------------
inline string ModelTree::getExpression(NodeID StartID, EquationType  iEquationType, int iEquationID)
{

  // Stack of tokens
  stack <int, vector<NodeID> > stack_token;
  ostringstream exp;
  NodeID current_token_ID;

  stack_token.push(StartID);
  int precedence_last_op = 0;
  int last_op_code = 0;
  int on_the_right_of_upper_node = 0;

  while (stack_token.size() > 0)
    {
      current_token_ID = stack_token.top();
      // defining short hand for current op code
      int current_op_code = current_token_ID->op_code;
      if ( current_token_ID->tmp_status == 1)
        {
          exp << "T" << current_token_ID->idx;
          // set precedence of terminal token to highest value
          precedence_last_op = 100;
          stack_token.pop();
          continue;
        }
      // else if current token is final
      else if ( current_op_code == NoOpCode)
        {
          exp << getArgument(current_token_ID->id1, current_token_ID->type1, iEquationType);
          // set precedence of terminal token to highest value
          precedence_last_op = 100;
          stack_token.pop();
          continue;
        }

      int precedence_current_op = operator_table.precedence(current_op_code);
      // deal with left argument first
      if (current_token_ID->left_done == 0)
        {
          // if current operator is not a temporary variable and
          // of lesser precedence than previous one, insert '('
          if ( precedence_current_op < precedence_last_op ||
               (on_the_right_of_upper_node == 1 &&
                (last_op_code == MINUS || last_op_code == DIVIDE) &&
                (precedence_current_op == precedence_last_op))||
               current_op_code == UMINUS)
            {
              exp << "(";
              current_token_ID->close_parenthesis = 1;
            }
          // set flag: left argument has been explored
          current_token_ID->left_done = 1;
          precedence_last_op = precedence_current_op;
          last_op_code = current_op_code;
          if ( offset == 0 && current_op_code == POWER)
            {
              exp << "pow(";
              precedence_last_op = 0;
              current_token_ID->close_parenthesis = 1;
            }
          else if ( current_op_code == UMINUS)
            {
              exp << "-";
              current_token_ID->close_parenthesis = 1;
            }
          else if ( operator_table.isfunction(current_op_code) == true)
            {
              exp << current_token_ID->op_name << "(";
              precedence_last_op = 0;
              current_token_ID->close_parenthesis = 1;
            }
          on_the_right_of_upper_node = 0;
          stack_token.push(current_token_ID->id1);
        }
      // deal with right argument when left branch is entirely explored
      else if ( current_token_ID->right_done == 0 )
        {
          current_token_ID->right_done = 1;
          if ( offset == 0 && current_op_code == POWER)
            {
              exp << ",";
            }

          if ( current_token_ID->id2 != NullID )
            {
              exp << current_token_ID->op_name;
              precedence_last_op = precedence_current_op;
              last_op_code = current_op_code;
              on_the_right_of_upper_node = 1;
              stack_token.push(current_token_ID->id2);
            }
        }
      else
        {
          if ( current_token_ID->close_parenthesis == 1)
            {
              exp << ")";
            }
          precedence_last_op = precedence_current_op;
          last_op_code = current_op_code;
          current_token_ID->left_done=0;
          current_token_ID->right_done=0;
          current_token_ID->close_parenthesis=0;
          stack_token.pop();
        }
    }
  return exp.str();
}

//------------------------------------------------------------------------------
inline string ModelTree::getArgument(NodeID id, Type type, EquationType iEquationType)
{

  stringstream  argument;

  if (type == eParameter)
    {
      argument << param_name << lpar << (long int)id+offset << rpar;
    }
  else if (type == eLocalParameter)
    {
      argument << SymbolTable::getNameByID(eLocalParameter,(long int)id);
    }
  else if (type == eNumericalConstant)
    {
      argument << NumericalConstants::get((long int) id);
    }
  else if (type == eEndogenous || type == eExogenous || type == eExogenousDet)
    if (iEquationType == eStaticEquations || iEquationType == eStaticDerivatives)
      {
        int idx = VariableTable::getSymbolID((long int) id)+offset;
        if (type == eEndogenous)
          {
            argument <<  "y" << lpar << idx << rpar;
          }
        else if (type == eExogenous)
          {
            argument << "x" << lpar << idx << rpar;
          }
        else if (type == eExogenousDet)
          {
            idx += ModelParameters::exo_nbr;
            argument <<  "x" << lpar << idx << rpar;
          }
      }
    else
      {
        if (type == eEndogenous)
          {
            int idx = VariableTable::getPrintIndex((long int) id)+offset;
            argument <<  "y" << lpar << idx << rpar;
          }
        else if (type == eExogenous)
          {
            int idx = VariableTable::getSymbolID((long int) id)+offset;
            int lag = VariableTable::getLag((long int) id);
            if (offset == 1)
              {
                if ( lag != 0)
                  {
                    argument <<  "x" << lpar << "it_ + " << lag
                             << ", " << idx << rpar;
                  }
                else
                  {
                    argument <<  "x" << lpar << "it_, " << idx << rpar;
                  }
              }
            else
              {
                if ( lag != 0)
                  {
                    argument <<  "x" << lpar << "it_+" << lag
                             << "+" << idx << "*nb_row_x" << rpar;
                  }
                else
                  {
                    argument <<  "x" << lpar << "it_+" << idx << "*nb_row_x" << rpar;
                  }
              }
          }
        else if (type == eExogenousDet)
          {
            int idx = VariableTable::getSymbolID((long int) id)+ModelParameters::exo_nbr+offset;
            int lag = VariableTable::getLag((long int) id);
            if (offset == 1)
              {
                if (lag != 0)
                  {
                    argument <<  "x" << lpar << "it_ + " << lag
                             << ", " << idx << rpar;
                  }
                else
                  {
                    argument <<  "x" << lpar << "it_, " << idx << rpar;
                  }
              }
            else
              {
                if (lag != 0)
                  {
                    argument <<  "x" << lpar << "it_ + " << lag
                             << "+" << idx <<  "*nb_row_xd" << rpar;
                  }
                else
                  {
                    argument <<  "x" << lpar << "it_+" << idx << "*nb_row_xd" <<  rpar;
                  }
              }

          }
      }
  return argument.str();
}

//------------------------------------------------------------------------------
void ModelTree::ModelInitialization(void)
{
  // Exit if there is no equation in model file*/
  if (ModelParameters::eq_nbr == 0)
    {
      (* error) ("no equations found in model file");
    }
  cout << ModelParameters::eq_nbr << " equation(s) found \n";
  // Sorting variable table
  VariableTable::Sort();

  // Setting number of equations in ModelParameters class
  // Here no derivative are computed
  BeginModel++;
  min_cost = 40*operator_table.cost(PLUS,offset);
  // Setting format of parentheses
  if (offset == 1)
    {
      lpar = '(';
      rpar = ')';
      param_name = "params";
    }
  else
    {
      lpar = '[';
      rpar = ']';
      param_name = "params";
    }
  /* Writing initialisation for M_.lead_lag_incidence matrix
     M_.lead_lag_incidence is a matrix with as many columns as there are
     endogenous variables and as many rows as there are periods in the
     models (nbr of rows = M_.max_lag+M_.max_lead+1)

     The matrix elements are equal to zero if a variable isn't present in the
     model at a given period.
  */
  // Initializing matrix to zero
  output << "M_.lead_lag_incidence = [";
  /*
    zeros(" <<
    ModelParameters::max_lag+ModelParameters::max_lead+1 << ", " <<
    ModelParameters::endo_nbr << ");\n";
  */
  // Loop on endogenous variables
  for (int endoID = 0; endoID < ModelParameters::endo_nbr; endoID++)
    {
      output << "\n\t";
      // Loop on periods
      for (int lag = -ModelParameters::max_endo_lag; lag <= ModelParameters::max_endo_lead; lag++)
        {
          // Getting name of symbol
          string name = SymbolTable::getNameByID(eEndogenous, endoID);
          // and its variableID if exists with current period
          int varID = VariableTable::getID(name, lag);
          //cout << name << " " << varID << " " << lag << " " << VariableTable::getPrintIndex(varID)+1 << " " << VariableTable::getSortID(varID)+1 << endl;

          if (varID >=0)
            {
              output << " " << VariableTable::getPrintIndex(varID)+1;
            }
          else
            {
              output << " 0";
            }
        }
      output << ";";
    }
  output << "]';\n";

  // Writing initialization for some other variables
  output << "M_.exo_names_orig_ord = [1:" << ModelParameters::exo_nbr << "];\n";
  output << "M_.maximum_lag = " << ModelParameters::max_lag << ";\n";
  output << "M_.maximum_lead = " << ModelParameters::max_lead<< ";\n";
  if (ModelParameters::endo_nbr)
    {
      output << "M_.maximum_endo_lag = " << ModelParameters::max_endo_lag << ";\n";
      output << "M_.maximum_endo_lead = " << ModelParameters::max_endo_lead<< ";\n";
      output << "oo_.steady_state = zeros(" << ModelParameters::endo_nbr << ", 1);\n";
    }
  if (ModelParameters::exo_nbr)
    {
      output << "M_.maximum_exo_lag = " << ModelParameters::max_exo_lag << ";\n";
      output << "M_.maximum_exo_lead = " << ModelParameters::max_exo_lead<< ";\n";
      output << "oo_.exo_steady_state = zeros(" << ModelParameters::exo_nbr << ", 1);\n";
    }
  if (ModelParameters::exo_det_nbr)
    {
      output << "M_.maximum_exo_det_lag = " << ModelParameters::max_exo_det_lag << ";\n";
      output << "M_.maximum_exo_det_lead = " << ModelParameters::max_exo_det_lead<< ";\n";
      output << "oo_.exo_det_steady_state = zeros(" << ModelParameters::exo_det_nbr << ", 1);\n";
    }
  if (ModelParameters::recur_nbr)
    {
      output << "M_.maximum_recur_lag = " << ModelParameters::max_recur_lag << ";\n";
      output << "M_.maximum_recur_lead = " << ModelParameters::max_recur_lead<< ";\n";
      output << "oo_.recur_steady_state = zeros(" << ModelParameters::recur_nbr << ", 1);\n";
    }
  if (ModelParameters::parameter_nbr)
    {
      output << "M_.params = zeros(" << ModelParameters::parameter_nbr << ", 1);\n";
    }
}

//------------------------------------------------------------------------------
string ModelTree::get()
{
  return output.str();
}

//------------------------------------------------------------------------------
inline int ModelTree::optimize(NodeID node)
{
  int cost;
  int tmp_status = 0;
  if (node->op_code != NoOpCode)
    {
      cost = operator_table.cost(node->op_code,offset);
      if (node->id1 != NullID && node->id1->op_code != NoOpCode)
        {
          cost += operator_table.cost(node->id1->op_code,offset);
        }
      if (node->id2 != NullID && node->id2->op_code != NoOpCode)
        {
          cost += operator_table.cost(node->id2->op_code,offset);
        }
    }
  cost *= node->reference_count[current_order];
  if (cost > min_cost)
    {
      tmp_status = 1;
      node->cost = 0;
    }
  else
    {
      tmp_status = 0;
      node->cost = cost;
    }
  node->tmp_status = 0;
  return tmp_status;
}

//------------------------------------------------------------------------------
#ifdef TEST_MODELTREE
int main(void)
{
  SymbolTable     st;
  VariableTable     vt;
  NumericalConstants  nc;
  ModelTree     model;
  vector<int>     t(20);

  //Adding 2 different symbols with AddSymbolDeclar
  SymbolTable::AddSymbolDeclar("c",eExogenous);
  //SymbolTable::PrintSymbolTable();
  SymbolTable::AddSymbolDeclar("k",eEndogenous);
  //SymbolTable::PrintSymbolTable();
  SymbolTable::AddSymbolDeclar("aa",eParameter);
  //SymbolTable::PrintSymbolTable();
  SymbolTable::AddSymbolDeclar("x",eExogenous);
  SymbolTable::AddSymbolDeclar("alph",eParameter);
  SymbolTable::AddSymbolDeclar("delt",eParameter);

  VariableTable::AddVariable("k",0);
  VariableTable::AddVariable("x",-1);
  VariableTable::AddVariable("c",-1);

  SymbolTable::AddSymbolDeclar("x1",eEndogenous);
  SymbolTable::AddSymbolDeclar("x2",eExogenousDet);
  //SymbolTable::AddSymbolDeclar("x3",eExogenous);

  VariableTable::AddVariable("x1",-1);
  VariableTable::AddVariable("x2",1);
  //VariableTable::AddVariable("x3",-1);
  //VariableTable::AddVariable("k",1);
  //VariableTable::AddVariable("y",0);

  t[0] = model.AddToken("aa");
  t[1] = model.AddToken("x",-1);
  t[2] = model.AddToken("k",0);
  t[3] = model.AddToken(Argument(t[0], eTempResult),
                        Argument(t[1], eTempResult), TIMES);
  t[4] = model.AddToken("alph");
  t[5] = model.AddToken(Argument(t[2], eTempResult),
                        Argument(t[4], eTempResult), POWER);
  t[6] = model.AddToken(Argument(t[3], eTempResult),
                        Argument(t[5], eTempResult), TIMES);

  t[7] = model.AddToken("delt");
  //t[8] = model.AddToken("1");
  t[9] = model.AddToken(Argument(1, eTempResult),
                        Argument(t[7], eTempResult), MINUS);
  t[10] = model.AddToken(Argument(t[9], eTempResult),
                         Argument(t[2], eTempResult), TIMES);

  t[11] = model.AddToken(Argument(t[2], eTempResult),
                         UMINUS);

  t[12] = model.AddToken(Argument(t[11], eTempResult),
                         Argument(t[6], eTempResult), PLUS);
  t[13] = model.AddToken(Argument(t[12], eTempResult),
                         Argument(t[10], eTempResult), PLUS);
  t[14] = model.AddToken(Argument(t[13], eTempResult),
                         Argument(t[10], eTempResult), PLUS);
  t[15] = model.AddToken("c",-1);
  t[16] = model.AddToken(Argument(t[15], eTempResult),
                         Argument(t[14], eTempResult), EQUAL);
  //try
  //{
  model.derive(2);
  model.setStaticModel();
  model.setDynamicStochasticModel();
  model.Open("static_model.m", "dynamic_model.m");
  model.Save();
  //cout << model.getStaticModel();
  //}
  //catch(Error err)
  //{
  // cout << "error---------------------\n";
  // exit(-1);
  //}
  //cout << model.getDynamicDeterministicModel() << endl;
  //cout << model.getDynamicStochasticModel() << endl;
  //VariableTable::Sort();
}
#endif
//------------------------------------------------------------------------------
