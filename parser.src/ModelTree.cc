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
#include "VariableTable.h"
#include "d_tab.h"
#include "NumericalConstants.h"
#include "ModelTree.h"
//------------------------------------------------------------------------------
ostringstream ModelTree::output;
//------------------------------------------------------------------------------
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
		iModelFileName1 += ".m"; 
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
		mStaticModelFile << "%\n% Status : Computes static model for Dynare\n%\n";
		mStaticModelFile << "% Warning : this file is generated automatically by Dynare\n";
		mStaticModelFile << "%			from model file (.mod)\n\n";
		if (iModelFileName2.size() && (computeJacobian||computeJacobianExo||computeHessian))
		{
			iModelFileName2 += ".m"; 
			mDynamicModelFile.open(iModelFileName2.c_str(),ios::out|ios::binary);
			if (!mDynamicModelFile.is_open())
			{
				cout << "ModelTree::Open : Error : Can't open file " << iModelFileName2
					<< " for writing\n";
				exit(-1);
			}
			iModelFileName2.erase(iModelFileName2.end()-2,iModelFileName2.end());
			mDynamicModelFile << "function [residual, g1, g2] = " <<  iModelFileName2 << "(y, x)\n";
			mDynamicModelFile << "%\n% Status : Computes dynamic model for Dynare\n%\n";
			mDynamicModelFile << "%Warning : this file is generated automatically by Dynare\n";
			mDynamicModelFile << "%			from model file (.mod)\n\n";
		
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
		mStaticModelFile << " *	          from model file (.mod)\n\n";
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
			mDynamicModelFile << " *	          from model file (.mod)\n\n";
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
		mStaticModelFile.close();
	}
	if (mDynamicModelFile.is_open() && (computeJacobian||computeJacobianExo||computeHessian))
	{
		mDynamicModelFile << DynamicOutput.str();
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
	NodeID	lToken;				// To store current working token
	NodeID	 	lD1, lD2;			// To store derivative arguments of 
									// current argument
	NodeID		lArg1, lArg2;		// To store	current arguments
	Type 		lType1;				// Type of first argument
	NodeID		t1,t11,t12,t13,
				t14, t15; 			// To store temoporary result arguments
	TreeIterator	BeginIT;		// Iterator of the 1st token to derive
	TreeIterator	EndIT;			// Iterator of the last token to derive
	TreeIterator	currentIT;		// Iterator (counter) for model tree loop

	vector<NodeID> EqualTokenIDs;		// IDs of "equal token" in model Tree
	// Capturing equation IDs
	for (currentIT = BeginModel; currentIT != mModelTree.end(); currentIT++)
	{	
		if ((*currentIT)->op_code == EQUAL)
		{
			EqualTokenIDs.push_back(*currentIT);
			// Equation is forced to be in Model Tree as refferenced
			// This is usfull to remove symetric elements
			(*currentIT)->reference_count[0]++;
		}
	} 
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
	// initialize derivatives of variables 
	int nbr_deriv_var = 0;

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
		// Loop on variables of derivation
		for (int var = 0; var < VariableTable::size(); var++)
		{		
			
			
			// Loop on tokens			
	  		for (currentIT = BeginIT;; currentIT++)
			{	
			  //cout << "Token " << (*currentIT)->idx << endl;
			  if (accumulate((*currentIT)->reference_count.begin(), (*currentIT)->reference_count.end(),0) > 0)
			  {
				lToken = *currentIT;//mModelTree[TokenCount];

				lArg1 = lToken->id1;
				lArg2 = lToken->id2;
				lType1 = lToken->type1;
				lD1 = DeriveArgument(lArg1, lType1, var);
				if (lArg2 != NullID) 
					lD2 = DeriveArgument(lArg2, eTempResult, var);
				// Case where token is a final argument
				if (lToken->op_code == NoOpCode)
				{
						setDerivativeAdress(*currentIT, lD1, var);
				}
				else
				{

				  switch (lToken->op_code)
				  {
				  	case UMINUS:
				  		t1 = AddUMinus(lD1);
						setDerivativeAdress(*currentIT, t1, var);
					 	break;
					case PLUS:
						t1 = AddPlus(lD1, lD2);
						setDerivativeAdress(*currentIT, t1, var);
					 	break;
					case MINUS:
						t1 = AddMinus(lD1, lD2);
						setDerivativeAdress(*currentIT, t1, var);
					 	break;
					case TIMES:
						t11 = AddTimes(lD1, lArg2);
						t12 = AddTimes(lD2, lArg1);
						t1 = AddPlus(t11, t12);
						setDerivativeAdress(*currentIT, t1, var);
						break;
					case DIVIDE:
						t11 = AddTimes(lD1, lArg2);
						t12 = AddTimes(lD2, lArg1);
						t13 = AddMinus(t11, t12);
						t14 =  AddTimes(lArg2, lArg2);
						t1 = AddDivide(t13, t14);
						setDerivativeAdress(*currentIT, t1, var);
						break;	
					case SQRT:
						t11 = AddPlus(*currentIT, *currentIT);
						t1 = AddDivide(lD1, t11);
						setDerivativeAdress(*currentIT, t1,var);
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
						setDerivativeAdress(*currentIT, t1, var);
						break;
					case EXP:
						t1 = AddTimes(lD1, *currentIT);
						setDerivativeAdress(*currentIT, t1, var);
						break;
					case LOG:
						t1 = AddDivide(lD1, lArg1);
						setDerivativeAdress(*currentIT, t1, var);
						break;
					case LOG10:
						t11 = AddExp(One);
						t12 = AddLog10(t11);
						t13 = AddDivide(lD1, lArg1);
						t1 = AddTimes(t12, t13);
						setDerivativeAdress(*currentIT, t1, var);
						break;
					case COS:
						t11 = AddSin(lArg1);
						t12 = AddUMinus(t11);
						t1 = AddTimes( lD1, t12);
						setDerivativeAdress(*currentIT, t1, var);
						break;
					case SIN:
						t11 = AddCos(lArg1);
						t1 = AddTimes(lD1,t11);
						setDerivativeAdress(*currentIT, t1, var);
						break;
					case TAN:
						t11 = AddTimes(*currentIT, *currentIT);
						t12 = AddPlus(t11, One);
						t1 = AddTimes(lD1, t12);
						setDerivativeAdress(*currentIT, t1, var);
		      			break;
					case ACOS:
						t11 = AddSin(*currentIT);
						t12 = AddDivide(lD1, t11);
						t1 = AddUMinus(t12);
						setDerivativeAdress(*currentIT, t1, var);
						break;
					case ASIN:
						t11 = AddCos(*currentIT);
						t1 = AddDivide(lD1, t11);
						setDerivativeAdress(*currentIT, t1, var);
						break;
					case ATAN:
						t11 = AddTimes(lArg1, lArg1);
						t12 = AddPlus(One, t11);
						t1 = AddDivide(lD1, t12);
						setDerivativeAdress(*currentIT, t1,var);
						break;
					case COSH:
						t11 = AddSinH(lArg1);
						t1 = AddTimes( lD1,t11);
						setDerivativeAdress(*currentIT, t1, var);
						break;
					case SINH:
						t11 = AddCosH(lArg1);
						t1 = AddTimes( lD1, t11);
						setDerivativeAdress(*currentIT, t1, var);
						break;
					case TANH:
						t11 = AddTimes(*currentIT, *currentIT);
						t12 = AddMinus(One, t11);
						t1 = AddTimes(lD1, t12);
						setDerivativeAdress(*currentIT, t1,var);
						break;
					case ACOSH:
						t11 = AddSinH(*currentIT);
						t1 = AddDivide(lD1, t11);
						setDerivativeAdress(*currentIT, t1,var);
						break;
					case ASINH:
						t11 = AddCosH(*currentIT);
						t1 = AddDivide(lD1, t11);
						setDerivativeAdress(*currentIT, t1,var);
						break;
		      		case ATANH:
						t11 = AddTimes(lArg1, lArg1);
						t12 = AddMinus(One, t11);
						t1 = AddTimes(lD1, t12);
						setDerivativeAdress(*currentIT, t1,var);
						break;
				  }
				}
			  }
			  if (currentIT == EndIT) break;
			}

			// Treating equal tokens
			// Skeeping symetric elements
			//vector<MetaToken>::iterator tree_it2 = mModelTree.begin();
			//int id = 0;
			int starti = var*Order*(Order-1)*ModelParameters::eq_nbr/2;
	  		for (int i = starti; i < EqualTokenIDs.size() ; i++ )
			{
				lToken = EqualTokenIDs[i];
				lArg1 = lToken->id1;
				lArg2 = lToken->id2;
				lType1 = lToken->type1;
				lD1 = DeriveArgument(lArg1, lType1, var);
				lD2 = DeriveArgument(lArg2, eTempResult, var);
				// If one hand sid is null, take the other 
				if (lD1 == Zero && lD2 != Zero)
				{
					t11 = AddUMinus(lD2);
					t1 = AddEqual(t11, Zero);
				}
				else if (lD1 != Zero && lD2 == Zero)
				{
					t1 = AddEqual(lD1, Zero);
				}
				else
				{
					t11 = AddMinus(lD1, lD2);
					t1 = AddEqual(t11, Zero);
				}
				// The derivative is forced to be in Model Tree as refferenced
				// This is usfull to remove symetric elements
				IncrementReferenceCount(t1);
				setDerivativeAdress(EqualTokenIDs[i], t1, var);
				if (Order == 1)
				{
					mDerivativeIndex[0].push_back(DerivativeIndex(t1, i-starti, var));					
				}
				else if (Order == 2)
				{
					int var1 = VariableTable::getSortID(i/ModelParameters::eq_nbr);
					int var2 = VariableTable::getSortID(var);
					mDerivativeIndex[1].push_back(DerivativeIndex(
						t1,
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
		for (int i=0; i< mDerivativeIndex[Order-1].size();i++)
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
		cout << "done \n";

	}
}
//------------------------------------------------------------------------------
inline bool ModelTree::writeAsTemp(NodeID id)
{
	int ref_count = id->reference_count[current_order];
	if (( ref_count > 1) && (ref_count*(id->cost) > min_cost))
		return true;
	else
		return false;
}
//------------------------------------------------------------------------------
inline NodeID ModelTree::DeriveArgument(NodeID iArg, Type iType, int iVarID)
{
	NodeID d;
	switch(iType)
	{
		case eTempResult			:
				d = iArg->d1[iVarID];
				return d;
			break;
		case eExogenous 			:
		case eExogenousDet 			:
		case eEndogenous			:
		case eRecursiveVariable		:
			if ((int) iArg == iVarID)
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
		case eNumericalConstant  	:
		case eParameter 			: 
			return Zero; 
			break;
		case eUNDEF					:
			return NullID;
			break;
		case eLoopIndex			:
		case eUnkownFunction		:
			return Zero;
			break;
		default				:
			cout << "ModelTree::DeriveArgument : Error :Unkown Type \n";
			exit(-1);
	};
	
}
//------------------------------------------------------------------------------
inline void ModelTree::setDerivativeAdress(NodeID iTokenID, NodeID iDerivative,int iVarID)
{
	//Derivative	lDerivative;
	
	//lDerivative.variable_id = iVarID;
	//lDerivative.derivative_address = iDerivative;
	iTokenID->d1[iVarID] = iDerivative;
	//mModelTree[iDerivative].p1[iVarID] = iTokenID;
}
//------------------------------------------------------------------------------
string 	ModelTree::setStaticModelM(void)
{
	TreeIterator tree_it;
	int lEquationNBR = 0;
	ostringstream model_output;		// Used for storing model equations
	ostringstream model_tmp_output;	// Used for storing tmp expressions for model equations
	ostringstream jacobian_output;	// Used for storing jacobian equations
	ostringstream jacobian_tmp_output;	// Used for storing tmp expressions for jacobian equations
	
	int	d = current_order;  		// Minimum number of times a temparary expression apears in equations
	int EquationNBR;				// Number of model equations
	int col = 1;					// Colomn index of Jacobian
		
	EquationNBR = ModelParameters::eq_nbr;
	// Reference count of token "0=0" is set to 0
	// Not to be printed as a temp expression
  	fill(ZeroEqZero->reference_count.begin(),
					   ZeroEqZero->reference_count.end(),0);
	// Setting tmp_status to 0, 
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		(*tree_it)->tmp_status = 0;
		
	}
	// Writing model Equations
	current_order = 1;
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if ((*tree_it)->op_code == EQUAL)
		{
			if  (lEquationNBR < ModelParameters::eq_nbr) 
			{
				model_output << getExpression(*tree_it, eStaticEquations, lEquationNBR) << endl;;
				lEquationNBR++;
			}
			else break;
		}			
	}

	for (tree_it = BeginModel; tree_it != mModelTree.end();tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			model_tmp_output << "  T" << (*tree_it)->idx << " = " << (*tree_it)->exp[0] << ";\n";
		}
	}

	// Writing Jacobian for endogenous variables without lag
	lEquationNBR = 0;
	for (int i = 0; i < mDerivativeIndex[0].size(); i++)
	{
		if (VariableTable::getType(mDerivativeIndex[0][i].derivators) == eEndogenous)
		{			
			NodeID startJacobian = mDerivativeIndex[0][i].token_id;
			string exp = getExpression(startJacobian, eStaticDerivatives);
			if (startJacobian != ZeroEqZero)
			{
				ostringstream g1;
				g1 << "  g1(" << mDerivativeIndex[0][i].equation_id+1 << ", " <<
					 VariableTable::getSymbolID(mDerivativeIndex[0][i].derivators)+1 << ')'; 
				jacobian_output << g1.str() << "=" <<  g1.str() << "+" << exp << ";\n";
			}
		}
	}
	for (tree_it = BeginModel; tree_it != mModelTree.end();tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			jacobian_tmp_output << "  T" << (*tree_it)->idx << " = " << (*tree_it)->exp[0] << ";\n";
		}
	}

	// Writing ouputs
	StaticOutput << "global M_ \n";

	StaticOutput << "if nargout >= 1,\n";
	StaticOutput << "  residual = zeros( " << ModelParameters::eq_nbr << ", 1);\n";
	StaticOutput << "\n\t%\n\t% Model equations\n\t%\n\n";
    StaticOutput << model_tmp_output.str() << "\n";
    StaticOutput << model_output.str();
    StaticOutput << "end\n";
	StaticOutput << "if nargout >= 2,\n";
	StaticOutput << "  g1 = " << 
		"zeros(" << ModelParameters::eq_nbr << ", " <<
		ModelParameters::endo_nbr << ");\n" ;
	StaticOutput << "\n\t%\n\t% Jacobian matrix\n\t%\n\n";
    StaticOutput << jacobian_tmp_output.str() << "\n";
	StaticOutput << jacobian_output.str();		
	StaticOutput << "end\n";
	current_order = d;
	return StaticOutput.str();
}
//------------------------------------------------------------------------------
string  ModelTree::setDynamicModelM(void)
{
	TreeIterator tree_it;
	int lEquationNBR = 0;
	ostringstream lsymetric;			// Used when writing symetric elements
	ostringstream model_output;			// Used for storing model equations
	ostringstream model_tmp_output;		// Used for storing tmp expressions for model equations
	ostringstream jacobian_output;		// Used for storing jacobian equations
	ostringstream jacobian_tmp_output;	// Used for storing tmp expressions for jacobian equations
	ostringstream hessian_output;		// Used for storing Hessian equations
	ostringstream hessian_tmp_output;	// Used for storing tmp expressions for Hessian equations
	
	int d = current_order;
	
 
  
  // Reference count of token "0=0" is set to 0
  // Not to be printed as a temp expression
  fill(ZeroEqZero->reference_count.begin(),
					   ZeroEqZero->reference_count.end(),0);
  // Setting tmp_status to 0, 
  for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
  {
		(*tree_it)->tmp_status = 0;
		
  }

  DynamicOutput << "global M_ it_\n";

  // Case where residuals and Jacobian with respect to endogenous variables
  // are computed
  if (computeJacobian && !computeJacobianExo && !computeHessian)
  {
  	DynamicOutput << "g2 = " << 
		"zeros(" << ModelParameters::eq_nbr << ", " <<
		ModelParameters::var_endo_nbr*ModelParameters::var_endo_nbr << ");\n" ;

	// Clearing output string
	model_output.str("");
	jacobian_output.str("");
	model_tmp_output.str("");
	jacobian_tmp_output.str("");
	// Getting equations from model tree
	// Starting from the end of equation
	// Searching for the next '=' operator,
	current_order = 1;
	lEquationNBR = 0;
	cout << "\tequations .. ";
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if ((*tree_it)->op_code == EQUAL)
		{
			if  (lEquationNBR < ModelParameters::eq_nbr) 
			{
				model_output << getExpression(*tree_it, eDynamicEquations, lEquationNBR) << endl;;
				lEquationNBR++;
			}
			else break;
		}			
	}
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			model_tmp_output << "  T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}

	cout << "done \n";
	// Getting Jacobian from model tree
	cout << "\tJacobian .. ";
	lEquationNBR = 0;
	for (int i = 0; i < mDerivativeIndex[0].size(); i++)
	{
		if (VariableTable::getType(mDerivativeIndex[0][i].derivators) == eEndogenous)
		{			
			NodeID startJacobian = mDerivativeIndex[0][i].token_id;
			string exp = getExpression(startJacobian, eDynamicDerivatives);
			if (startJacobian != ZeroEqZero)
				jacobian_output << "  g1(" << mDerivativeIndex[0][i].equation_id+1 << ", " <<
					VariableTable::getPrintIndex(mDerivativeIndex[0][i].derivators)+1 << ") = " << exp << ";\n";
		}
	}
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			jacobian_tmp_output << "  T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}
	cout << "done \n";

	DynamicOutput << "\n%\n% Computting residuals and Jacobian with respect\n";
	DynamicOutput << "% to endogenous variables\n%\n";
	DynamicOutput << "if nargout >= 1,\n";
	DynamicOutput << "  residual = zeros( " << ModelParameters::eq_nbr << ", 1);\n";
	DynamicOutput << "\n\t%\n\t% Model equations\n\t%\n\n";
	DynamicOutput << model_tmp_output.str() << "\n";
	DynamicOutput << model_output.str();
	DynamicOutput << "end\n";
	DynamicOutput << "if nargout >= 2,\n";
	DynamicOutput << "\n\t%\n\t% Jacobian matrix\n\t%\n\n";
  	DynamicOutput << "  g1 = " << 
		"zeros(" << ModelParameters::eq_nbr << ", " <<
		ModelParameters::var_endo_nbr << ");\n" ;

	DynamicOutput << jacobian_output.str();
	DynamicOutput << jacobian_tmp_output.str() << "\n";
	DynamicOutput << "end\n";
  }
  // Case where residuals and Jacobian with respect to endogenous and exogenous 
  // variables are computed
  else if (computeJacobianExo && !computeHessian)
  {
	// Clearing output string
	model_output.str("");
	jacobian_output.str("");
	model_tmp_output.str("");
	jacobian_tmp_output.str("");
	current_order = 1;
	lEquationNBR = 0;
	// Getting equations from model tree
	// Starting from the end of equation
	// Searching for the next '=' operator,
	
	cout << "\tequations .. ";
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if ((*tree_it)->op_code == EQUAL)
		{
			if  (lEquationNBR < ModelParameters::eq_nbr) 
			{
				model_output << getExpression(*tree_it, eDynamicEquations, lEquationNBR) << endl;
				lEquationNBR++;           
			}                         
			else break;
		}			
	}
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			model_tmp_output << "  T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}
	cout << "done \n";

	// Getting Jacobian from model tree
	// Starting from the end of equation
	// Searching for the next '=' operator,
	
	cout << "\tJacobian .. ";
	for (int i = 0; i < mDerivativeIndex[0].size(); i++)
	{
		NodeID startJacobian = mDerivativeIndex[0][i].token_id;
		string exp = getExpression(startJacobian, eDynamicDerivatives);
		if (startJacobian != ZeroEqZero)
			jacobian_output << "  g1(" << mDerivativeIndex[0][i].equation_id+1 << ", " <<
				VariableTable::getSortID(mDerivativeIndex[0][i].derivators)+1 << ") = " << exp << ";\n";
	}
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			jacobian_tmp_output << "  T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}
	cout << "done \n";

	DynamicOutput << "\n%\n% Computting residuals and Jacobian with respect\n";
	DynamicOutput << "% to endogenous and exogenous variables\n%\n";
	DynamicOutput << "if nargout >= 1,\n";
	DynamicOutput << "  residual = zeros( " << ModelParameters::eq_nbr << ", 1);\n";
	DynamicOutput << "\n\t%\n\t% Model equations\n\t%\n\n";
    DynamicOutput << model_tmp_output.str() << "\n";
    DynamicOutput << model_output.str();
    DynamicOutput << "end\n";
    DynamicOutput << "if nargout >= 2,\n";
 	// Writing initialization instruction for matrix g1
	DynamicOutput << "  g1 = " << 
		"zeros(" << ModelParameters::eq_nbr << ", " <<
		VariableTable::size() << ");\n" ;
	DynamicOutput << "\n\t%\n\t% Jacobian matrix\n\t%\n\n";
	DynamicOutput << jacobian_tmp_output.str() << "\n";
	DynamicOutput << jacobian_output.str();
	DynamicOutput << "end\n";
  }
  // Case where residuals, Jacobian and Hessian with respect to endogenous and exogenous
  // variables are computed
  else if (computeHessian && computeJacobianExo)
  {
   	// Clearing output string
	model_output.str("");
	jacobian_output.str("");
	hessian_output.str("");
	model_tmp_output.str("");
	jacobian_tmp_output.str("");
	hessian_tmp_output.str("");
	lsymetric.str("");

	current_order = 2;
	lEquationNBR = 0;
	// Getting equations from model tree
	// Starting from the end of equation
	// Searching for the next '=' operator,
	
	cout << "\tequations .. ";
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if ((*tree_it)->op_code == EQUAL)
		{
			if  (lEquationNBR < ModelParameters::eq_nbr) 
			{
				model_output << getExpression(*tree_it, eDynamicEquations, lEquationNBR) << endl;
				lEquationNBR++;           
			}                         
			else break;
		}			
	}
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			model_tmp_output << "  T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}
	cout << "done \n";

	// Getting Jacobian from model tree
	// Starting from the end of equation
	// Searching for the next '=' operator,
	
	cout << "\tJacobian .. ";
	for (int i = 0; i < mDerivativeIndex[0].size(); i++)
	{
		NodeID startJacobian = mDerivativeIndex[0][i].token_id;
		string exp = getExpression(startJacobian, eDynamicDerivatives);
		if (startJacobian != ZeroEqZero)
			jacobian_output << "  g1(" << mDerivativeIndex[0][i].equation_id+1 << ", " <<
				VariableTable::getSortID(mDerivativeIndex[0][i].derivators)+1 << ") = " << exp << ";\n";
	}	
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			jacobian_tmp_output << "  T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}
	cout << "done \n";

	// Getting Hessian from model tree
	// Starting from the end of equation
	// Searching for the next '=' operator,		
 	lEquationNBR = 0;
	cout << "\tHessian .. ";
	for (int i = 0; i < mDerivativeIndex[1].size(); i++)
	{
		NodeID startHessian = mDerivativeIndex[1][i].token_id;
		string exp = getExpression(startHessian, eDynamicDerivatives);

		int varID1 = mDerivativeIndex[1][i].derivators/VariableTable::size();
		int varID2 = mDerivativeIndex[1][i].derivators-varID1*VariableTable::size();
		//cout << "ID = " << startHessian << " exp = " << exp << "\n";
		if (startHessian != ZeroEqZero)
		{
			hessian_output << "  g2(" << mDerivativeIndex[1][i].equation_id+1 << ", " <<
			  mDerivativeIndex[1][i].derivators+1 << ") = " << exp << ";\n";
			// Treating symetric elements
			if (varID1 != varID2)
	 			lsymetric <<  "  g2(" << mDerivativeIndex[1][i].equation_id+1 << ", " <<
						varID2*VariableTable::size()+varID1+1 << ") = " <<
						"g2(" << mDerivativeIndex[1][i].equation_id+1 << ", " <<
						mDerivativeIndex[1][i].derivators+1 << ");\n";
		}

	}
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			hessian_tmp_output << "  T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}
	cout << "done \n";

	DynamicOutput << "\n% Computting residuals, Jacobian and Hessian with respect\n";
	DynamicOutput << "% to endogenous and exogenous variables are computed\n";
	DynamicOutput << "if nargout >= 1,\n";
	DynamicOutput << "  residual = zeros( " << ModelParameters::eq_nbr << ", 1);\n";
	DynamicOutput << "\n\t%\n\t% Model equations\n\t%\n\n";
    DynamicOutput << model_tmp_output.str() << "\n";
    DynamicOutput << model_output.str();
    DynamicOutput << "end\n";
    DynamicOutput << "if nargout >= 2,\n";
	// Writing initialization instruction for matrix g2
	DynamicOutput << "  g1 = " << 
		"zeros(" << ModelParameters::eq_nbr << ", " <<
		VariableTable::size() << ");\n" ;
	DynamicOutput << "\n\t%\n\t% Jacobian matrix\n\t%\n\n";
	DynamicOutput << jacobian_tmp_output.str() << "\n";
	DynamicOutput << jacobian_output.str();
	DynamicOutput << "end\n";
	DynamicOutput << "if nargout >= 3,\n";
	// Writing initialization instruction for matrix g2
	DynamicOutput << "  g2 = " << 
		"zeros(" << ModelParameters::eq_nbr << ", " <<
		 VariableTable::size()*VariableTable::size()<< ");\n";
	DynamicOutput << "\n\t%\n\t% Hessian matrix\n\t%\n\n";
    DynamicOutput << hessian_tmp_output.str() << "\n";
    DynamicOutput << hessian_output.str() << lsymetric.str();
    DynamicOutput << "end;\n";
  }
  current_order = d;
  return DynamicOutput.str();
}
//------------------------------------------------------------------------------
string 	ModelTree::setStaticModelC(void)
{
	TreeIterator tree_it;
	int lEquationNBR = 0;
	ostringstream model_output;;		// Used for storing model equations
	ostringstream jacobian_output;		// Used for storing jacobian equations
	ostringstream model_tmp_output;;	// Used for storing tmp expressions formodel equations
	ostringstream jacobian_tmp_output;	// Used for storing tmp expressions for jacobian equations
	int	d = current_order;  			// Minimum number of times a temparary expression apears in equations
	int EquationNBR;					// Number of model equations
	int col = 1;						// Colomn index of Jacobian
		
	EquationNBR = ModelParameters::eq_nbr;
	// Reference count of token "0=0" is set to 0
	// Not to be printed as a temp expression
  	fill(ZeroEqZero->reference_count.begin(),
					   ZeroEqZero->reference_count.end(),0);
	// Setting tmp_status to 0, 
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		(*tree_it)->tmp_status = 0;
		
	}
	// Writing model Equations
	current_order = 1;
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if ((*tree_it)->op_code == EQUAL)
		{
			if  (lEquationNBR < ModelParameters::eq_nbr) 
			{
				model_output << getExpression(*tree_it, eStaticEquations, lEquationNBR) << endl;;
				lEquationNBR++;
			}
			else break;
		}			
	}
	for (tree_it = BeginModel; tree_it != mModelTree.end();tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			model_tmp_output << "  double T" << (*tree_it)->idx << " = " << (*tree_it)->exp[0] << ";\n";
		}
	}
	// Writing Jacobian for endogenous variables without lag
	lEquationNBR = 0;
	for (int i = 0; i < mDerivativeIndex[0].size(); i++)
	{
		if (VariableTable::getType(mDerivativeIndex[0][i].derivators) == eEndogenous)
		{			
			NodeID startJacobian = mDerivativeIndex[0][i].token_id;
			string exp = getExpression(startJacobian, eStaticDerivatives);
			if (startJacobian != ZeroEqZero)
			{
				ostringstream g1;
				g1 << "  g1[" << mDerivativeIndex[0][i].equation_id+
					 (VariableTable::getSymbolID(mDerivativeIndex[0][i].derivators))*ModelParameters::eq_nbr << "]"; 
				jacobian_output << g1.str() << "=" <<  g1.str() << "+" << exp << ";\n";
			}
		}
	}
	for (tree_it = BeginModel; tree_it != mModelTree.end();tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			jacobian_tmp_output << "  double T" << (*tree_it)->idx << " = " << (*tree_it)->exp[0] << ";\n";
		}
	}
	
	// Writing C function Static
	StaticOutput << "void Static(double *y, double *x, double *residual, double *g1)\n";
	StaticOutput << "{\n";
	StaticOutput << "  double lhs, rhs;\n\n";
	// Writing residual equations
	StaticOutput << "  /* Residual equations */\n";
	StaticOutput << "  if (residual == NULL) return;\n";
	StaticOutput << " {\n";
	StaticOutput << model_tmp_output.str() << "\n";
	StaticOutput << model_output.str();
	// Writing Jacobian
	StaticOutput << "   /* Jacobian for endogenous variables without lag */\n";
	StaticOutput << "   if (g1 == NULL) return;\n";
	StaticOutput << " {\n";
	StaticOutput << jacobian_tmp_output.str() << "\n";
	StaticOutput << jacobian_output.str();
	StaticOutput << "  }\n";		
	StaticOutput << " }\n";		
	StaticOutput << "}\n\n";
	current_order = d;
	return StaticOutput.str();
}
//------------------------------------------------------------------------------
string  ModelTree::setDynamicModelC(void)
{
	TreeIterator tree_it;
	int lEquationNBR = 0;
	ostringstream lsymetric;			// Used when writing symetric elements
	ostringstream model_output;;		// Used for storing model equations
	ostringstream jacobian_output;		// Used for storing jacobian equations
	ostringstream model_tmp_output;;	// Used for storing tmp expressions for model equations
	ostringstream jacobian_tmp_output;	// Used for storing tmp expressions for jacobian equations
	ostringstream hessian_output;		// Used for storing Hessian equations
	ostringstream hessian_tmp_output;	// Used for storing tmp expressions for Hessian equations
	int d = current_order;
	
 
  
  // Reference count of token "0=0" is set to 0
  // Not to be printed as a temp expression
  fill(ZeroEqZero->reference_count.begin(),
					   ZeroEqZero->reference_count.end(),0);
  // Setting tmp_status to 0, 
  for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
  {
		(*tree_it)->tmp_status = 0;
		
  }
  // Case where residuals and Jacobian with respect to endogenous variables
  // are computed
  if (computeJacobian && !computeJacobianExo && !computeHessian)
  {
	// Clearing output string
	model_output.str("");
	// Getting equations from model tree
	// Starting from the end of equation
	// Searching for the next '=' operator,
	current_order = 1;
	lEquationNBR = 0;
	cout << "\tequations .. ";
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if ((*tree_it)->op_code == EQUAL)
		{
			if  (lEquationNBR < ModelParameters::eq_nbr) 
			{
				model_output << getExpression(*tree_it, eDynamicEquations, lEquationNBR) << endl;;
				lEquationNBR++;
			}
			else break;
		}			
	}
	for (tree_it = BeginModel;tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			model_tmp_output << " double T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}

	cout << "done \n";
	// Getting Jacobian from model tree
	cout << "\tJacobian .. ";
	lEquationNBR = 0;
	for (int i = 0; i < mDerivativeIndex[0].size(); i++)
	{
		if (VariableTable::getType(mDerivativeIndex[0][i].derivators) == eEndogenous)
		{			
			NodeID startJacobian = mDerivativeIndex[0][i].token_id;
			string exp = getExpression(startJacobian, eDynamicDerivatives);
			if (startJacobian != ZeroEqZero)
				jacobian_output << "  g1[" << mDerivativeIndex[0][i].equation_id+
					(VariableTable::getPrintIndex(mDerivativeIndex[0][i].derivators))*
					ModelParameters::eq_nbr << "] = " << exp << ";\n";
		}
	}
	for (tree_it = BeginModel;tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			jacobian_tmp_output << " double T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}	
	cout << "done \n";

	// Writing C function Dynamic
	DynamicOutput << "void Dynamic(double *y, double *x, double *residual, double *g1, double *g2)\n";
	DynamicOutput << "{\n";
	DynamicOutput << "  double lhs, rhs;\n\n";
	DynamicOutput << "  /* Residual equations */\n";
	DynamicOutput << "  if (residual == NULL) return;\n";
	DynamicOutput << " {\n";
	DynamicOutput << model_output.str() << "\n";
	DynamicOutput << model_output.str();
	DynamicOutput << "  /* Jacobian for for endogenous variables */\n";
	DynamicOutput << "  if (g1 == NULL) return;\n";
	DynamicOutput << "  {\n";	
	DynamicOutput << jacobian_tmp_output.str() << "\n";
	DynamicOutput << jacobian_output.str();
	DynamicOutput << "  }\n";
	DynamicOutput << " }\n";
	DynamicOutput << "}\n\n";
  }
  // Case where residuals and Jacobian with respect to endogenous and exogenous 
  // variables are computed
  else if (computeJacobianExo && !computeHessian)
  {
	// Clearing output string
	model_output.str("");
	jacobian_output.str("");
	current_order = 1;
	lEquationNBR = 0;
	// Getting equations from model tree
	// Starting from the end of equation
	// Searching for the next '=' operator,
	
	cout << "\tequations .. ";
	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if ((*tree_it)->op_code == EQUAL)
		{
			if  (lEquationNBR < ModelParameters::eq_nbr) 
			{
				model_output << getExpression(*tree_it, eDynamicEquations, lEquationNBR) << endl;
				lEquationNBR++;           
			}                         
			else break;
		}			
	}
	for (tree_it = BeginModel;tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			model_tmp_output << " double T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}	
	cout << "done \n";

	// Getting Jacobian from model tree
	// Starting from the end of equation
	// Searching for the next '=' operator,
	
	cout << "\tJacobian .. ";
	for (int i = 0; i < mDerivativeIndex[0].size(); i++)
	{
		NodeID startJacobian = mDerivativeIndex[0][i].token_id;
		string exp = getExpression(startJacobian, eDynamicDerivatives);
		if (startJacobian != ZeroEqZero)
			jacobian_output << "  g1[" << mDerivativeIndex[0][i].equation_id +
				(VariableTable::getSortID(mDerivativeIndex[0][i].derivators))*
				ModelParameters::eq_nbr << "] = " << exp << ";\n";
	}	
	for (tree_it = BeginModel;tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			jacobian_tmp_output << " double T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}	
	cout << "done \n";

	// Writing C function Dynamic
	DynamicOutput << "void Dynamic(double *y, double *x, double *residual, double *g1, double *g2)\n";
	DynamicOutput << "{\n";
	DynamicOutput << "  double lhs, rhs;\n\n";
	DynamicOutput << " /* Residual equations */\n";
	DynamicOutput << " if (residual == NULL) return;\n";
	DynamicOutput << " {\n";
    DynamicOutput << model_tmp_output.str() << "\n";
    DynamicOutput << model_output.str();
	DynamicOutput << "  /* Jacobian for endogenous and exogenous variables */\n";
	DynamicOutput << "  if (g1 == NULL) return;\n";
	DynamicOutput << "  {\n";	
	DynamicOutput << jacobian_tmp_output.str() << "\n";
	DynamicOutput << jacobian_output.str();
	DynamicOutput << "  }\n";
	DynamicOutput << " }\n";
	DynamicOutput << "}\n\n";
  }
  // Case where residuals, Jacobian and Hessian with respect to endogenous and exogenous
  // variables are computed
  else if (computeHessian && computeJacobianExo)
  {
   	// Clearing output string
	model_output.str("");
	jacobian_output.str("");
	hessian_output.str("");
	lsymetric.str("");

	current_order = 2;
	lEquationNBR = 0;
	// Getting equations from model tree
	// Starting from the end of equation
	// Searching for the next '=' operator,
	
	cout << "\tequations .. ";

	for (tree_it = BeginModel; tree_it != mModelTree.end(); tree_it++)
	{
		if ((*tree_it)->op_code == EQUAL)
		{
			if  (lEquationNBR < ModelParameters::eq_nbr) 
			{
				model_output << getExpression(*tree_it, eDynamicEquations, lEquationNBR) << endl;
				lEquationNBR++;           
			}                         
			else break;
		}			
	}
	for (tree_it = BeginModel;tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			model_tmp_output << " double T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}	
	cout << "done \n";

	// Getting Jacobian from model tree
	// Starting from the end of equation
	// Searching for the next '=' operator,
	
	cout << "\tJacobian .. ";
	for (int i = 0; i < mDerivativeIndex[0].size(); i++)
	{
		NodeID startJacobian = mDerivativeIndex[0][i].token_id;
		string exp = getExpression(startJacobian, eDynamicDerivatives);
		if (startJacobian != ZeroEqZero)
			jacobian_output << "  g1[" << mDerivativeIndex[0][i].equation_id+
				(VariableTable::getSortID(mDerivativeIndex[0][i].derivators))
				*ModelParameters::eq_nbr << "] = " << exp << ";\n";
	}	
	for (tree_it = BeginModel;tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			jacobian_tmp_output << " double T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}	
	cout << "done \n";

	// Getting Hessian from model tree
	// Starting from the end of equation
	// Searching for the next '=' operator,		
 	lEquationNBR = 0;
	cout << "\tHessian .. ";
	for (int i = 0; i < mDerivativeIndex[1].size(); i++)
	{
		NodeID startHessian = mDerivativeIndex[1][i].token_id;
		string exp = getExpression(startHessian, eDynamicDerivatives);

		int varID1 = mDerivativeIndex[1][i].derivators/VariableTable::size();
		int varID2 = mDerivativeIndex[1][i].derivators-varID1*VariableTable::size();
		//cout << "ID = " << startHessian << " exp = " << exp << "\n";
		if (startHessian != ZeroEqZero)
		{
			hessian_output << "  g2[" << mDerivativeIndex[1][i].equation_id+
			  mDerivativeIndex[1][i].derivators*
			  ModelParameters::eq_nbr << "] = " << exp << ";\n";
			// Treating symetric elements
			if (varID1 != varID2)
	 			lsymetric <<  "  g2[" << mDerivativeIndex[1][i].equation_id+
						(varID2*VariableTable::size()+varID1)*
						ModelParameters::eq_nbr << "] = " <<
						"g2[" << mDerivativeIndex[1][i].equation_id+
						(mDerivativeIndex[1][i].derivators)*
						ModelParameters::eq_nbr << "];\n";
		}

	}
	for (tree_it = BeginModel;tree_it != mModelTree.end(); tree_it++)
	{
		if (((*tree_it)->tmp_status == 1) &&
			writeAsTemp(*tree_it))
		{
			(*tree_it)->tmp_status = -1;
			hessian_tmp_output << " double T" << (*tree_it)->idx << " = " << (*tree_it)->exp[1] << ";\n";
		}
	}	
	cout << "done \n";

	// Writing C function Dynamic
	DynamicOutput << "void Dynamic(double *y, double *x, double *residual, double *g1, double *g2)\n";
	DynamicOutput << "{\n";
	DynamicOutput << "  double lhs, rhs;\n\n";
	DynamicOutput << " /* Residual equations */\n";
	DynamicOutput << " if (residual == NULL) return;\n";
	DynamicOutput << " {\n";
    DynamicOutput << model_tmp_output.str() << "\n";
    DynamicOutput << model_output.str();
	DynamicOutput << "  /* Jacobian for endogenous and exogenous variables */\n";
	DynamicOutput << "  if (g1 == NULL) return;\n";
	DynamicOutput << "  {\n";	
	DynamicOutput << jacobian_tmp_output.str() << "\n";
	DynamicOutput << jacobian_output.str();
	DynamicOutput << "  /* Hessian for endogenous and exogenous variables */\n";
	DynamicOutput << "  if (g2 == NULL) return;\n";
	DynamicOutput << "   {\n";	
    DynamicOutput << hessian_tmp_output.str() << "\n";
    DynamicOutput << hessian_output.str() << lsymetric.str();
	DynamicOutput << "   }\n";
	DynamicOutput << "  }\n";
	DynamicOutput << " }\n";
	DynamicOutput << "}\n\n";
  }
  current_order = d;
  return DynamicOutput.str();
}
//------------------------------------------------------------------------------
inline string ModelTree::getExpression(NodeID StartID, EquationType  iEquationType, int iEquationID)
{

	// Stack of temporary tokens
	stack <int, vector<NodeID> > stack_token;
	// Stack of temporary expressions
	stack <int, vector<string> > stack_expression;
	// Temporary output 
	ostringstream exp;
	// temporary variables for saving arguments and name oparator
	string argument1, argument2, op_name;
	// Current token ID
	NodeID currentTokenID;
	int current_op, last_op;	

	// Initialization of "followed" flags
	StartID->followed1 = false;
	StartID->followed2 = false;
	// Setting equation type for exp field
	int eq_type;
	switch (iEquationType)
	{
		case eStaticEquations:
		case eStaticDerivatives:
			eq_type = 0;
			break;
		default :
			eq_type = 1;
	}


	stack_token.push(StartID);
	current_op = last_op = stack_token.top()->op_code; 
	currentTokenID = StartID;
	
	// Main loop : 
	// Repeat for last token from the stack
	// (1) 	if argument is temporary result, and not yet followed,
	//			set it as followed (flag) and push corresponding token 
	//			on the token stack
	// (2) argument followed, or final argument
	//		(2.1) if argument is followed 
	//			- set argument1 (or argument2) by last expression on
	//			expression tack
	//			- pop last expression from expression stack
	//		(2.2) if final argument
	//			  set argument1 (or argument2) by final argument
	// (3) set op_name by last token from the token stack
	// (3) pop last token from the token stack
	// (4) write temporary expression (using argument1, argument2 
	//		and op_name) and push it on the expression stack
	// (5)

	while (stack_token.size() > 0)
	{		
	  
	  currentTokenID = stack_token.top();

	  //currentTokenID = mIndexOfTokens[Key((MToken) stack_token.top())];
	  // Testing if expression has been writen as temp result
	  if ( currentTokenID->exp[eq_type].size() == 0 )
	  {
	  	//if (currentTokenID != 0 && currentTokenID != 3)
  		//	cout << "currentTokenID = " << currentTokenID << endl;
		// First argument is a temporary result,
		// pushing token on token stack and setting that argument to be followed
		if ((currentTokenID->followed1 == false) &&
			(currentTokenID->type1 == eTempResult))
		{
			currentTokenID->followed1 = true;
			// Initialization of "followed" flags
			currentTokenID->id1->followed1 = false;
			currentTokenID->id1->followed2 = false;
			stack_token.push(currentTokenID->id1);
		}
		// Second argument has id >=0 (temporary result),
		// pushing token on stack and setting that argument to be followed 
		else if ((currentTokenID->followed2 == false) && 
			(currentTokenID->id2 != NullID))
		{
			currentTokenID->followed2 = true;
			// Initialization of "followed" flags
			currentTokenID->id2->followed1 = false;
			currentTokenID->id2->followed2 = false;
			stack_token.push(currentTokenID->id2);
		}
		// Writing expression
		else
		{		
			//cout << "currentTokenID = " << currentTokenID << endl;	
			// Final token
			if (currentTokenID->op_code == NoOpCode)
			{
				argument1 = getArgument(currentTokenID->id1, currentTokenID->type1, iEquationType);
				//cout << "Test 1 : argument1 = " << argument1 << endl;
				current_op = last_op;
				op_name = currentTokenID->op_name;
				exp.str("");
				stack_token.pop();
				// Saving last operator for the followed argument
				if (stack_token.size() > 0)
				{
					last_op = stack_token.top()->op_code;
				}
				else
				{
					last_op = current_op;
				}
				exp << 	argument1;
				currentTokenID->tmp_status = 0;
			}
			// Testing if unary or binary token
			// Binary operator
			else if (currentTokenID->id2 != NullID)
			{								
				argument2 = stack_expression.top();	
				//cout << "Test 2 : argument2 = " << argument2 << endl;	
				stack_expression.pop();				
				argument1 = stack_expression.top();
				//cout << "Test 2 : argument1 = " << argument1 << endl;
				current_op = currentTokenID->op_code;
				stack_expression.pop();
				op_name = currentTokenID->op_name;
				exp.str("");
				stack_token.pop();
				// Saving last operator for the followed argument
				if (stack_token.size() > 0)
				{
					last_op = stack_token.top()->op_code;
				}
				else
				{
					last_op = current_op;
				}
				if (operator_table.precedence(current_op) < operator_table.precedence(last_op))
				{
					// Comma operator, no parentheses
					// parentheses are writing with function operator
					if (current_op == COMMA)
					{
						exp << argument1 << op_name << argument2 ;
					}								
					else if ((offset == 1) || (current_op != POWER))
					{
						exp << 	'(' << argument1 << op_name << argument2 << ')';
					}
					else
					{
						exp << 	'(' << "pow(" << argument1 << ", " << argument2 << "))";
					}
					currentTokenID->tmp_status = 1;
				}
				else
				{
					if (current_op == EQUAL)
					{
						//Writing model equations
						switch(iEquationType)
						{
						  case eDynamicDerivatives :
						  case eStaticDerivatives :
						  	{
							exp << argument1;
							NodeID id1 = currentTokenID->id1;
							currentTokenID->tmp_status = id1->tmp_status;
							//cout << "current_order = " << current_order << " : " <<
							//mModelTree[currentTokenID].reference_count[current_order] << " :	" << exp.str() << endl;
							}
							break;
						  case eDynamicEquations :
							exp << "  lhs = " << argument1 << ";\n";
							exp << "  rhs = " << argument2 << ";\n";
							exp << "  residual" << lpar << iEquationID+offset << rpar << " = lhs - rhs;\n";
							currentTokenID->tmp_status = 0;
							//cout << "current_order = " << current_order << " : " 
							//<< mModelTree[currentTokenID].reference_count[current_order] << " :	" << exp.str();
							break;
						  case eStaticEquations : 
							exp << "  lhs = " << argument1 << ";\n";
							exp << "  rhs = " << argument2 << ";\n";
							exp << "  residual" << lpar << iEquationID+offset << rpar << " = lhs - rhs;\n";
							currentTokenID->tmp_status = 0;
							break;
						}
					}
					else
					{	// Matlab format
						if (offset == 1)
							exp << 	argument1 << op_name << argument2;
						// C format
						else
							if (current_op == POWER)
								exp << 	"pow(" << argument1 << "," << argument2 << ')';
							// In C language --X is not allowed
							else if (last_op == MINUS)
								exp << 	'(' << argument1 << op_name <<  argument2 << ')';
							else 
								exp << 	argument1 << op_name << argument2;
							
						currentTokenID->tmp_status = 1;
					}
				}
			}
			// Unary operator
			else
			{
				argument1 = stack_expression.top();
				//cout << "Test 2 : argument1 = " << argument1 << endl;
				current_op = currentTokenID->op_code;
				stack_expression.pop();	
				op_name = stack_token.top()->op_name;
				exp.str("");
				stack_token.pop();
				// Saving last operator for the followed argument
				if (stack_token.size() > 0)
				{
					last_op = stack_token.top()->op_code;
				}
				else
				{
					last_op = current_op;
				}
				/*
				if (operator_table.precedence(current_op) < operator_table.precedence(last_op))
				{
					exp << 	'(' << op_name << argument1 << ')';
					currentTokenID->tmp_status = 1;
					// Exclude "-cte" from temprary expressions
					if (currentTokenID->op_code == UMINUS)
					{
						NodeID id1 = currentTokenID->id1;
						if (id1->type1 == eNumericalConstant)
							currentTokenID->tmp_status = 0;
					}
						
				}	
				else
				{	
				*/
					currentTokenID->tmp_status = 1;
					// Case of functions 
					if (operator_table.isfunction(current_op) == true)
					{
						exp << 	op_name << '(' << argument1 << ')';
					}
					else
					{
						exp << 	op_name << '(' << argument1 << ')';
						// Exclude "-cte" from temprary expressions
						NodeID id1 = currentTokenID->id1;
						if (id1->type1 != eTempResult)
							currentTokenID->tmp_status = 0;
					}
				//}	
			}
			if (currentTokenID->tmp_status &&
				writeAsTemp(currentTokenID))
			{
				currentTokenID->exp[eq_type] = exp.str();
			  	exp.str("");
			  	exp << "T" << currentTokenID->idx;
				stack_expression.push(exp.str());
			}			
			else
			{
				stack_expression.push(exp.str());
				currentTokenID->exp[eq_type] = exp.str();

			}
		}
	  }
	  // Expression is in data member exp[]
	  else
	  {
	  	if (currentTokenID->tmp_status && 
	  		   writeAsTemp(currentTokenID))
		{
	  		exp.str("");
			stack_token.pop();
	  		exp << "T" << currentTokenID->idx;
	  		stack_expression.push(exp.str());
	  	}
	  	else
	  	{
	  		stack_token.pop();  		
	  		stack_expression.push(currentTokenID->exp[eq_type]);
		}
		last_op = current_op;
	  }

	}
	return stack_expression.top();
}
//------------------------------------------------------------------------------
/*
void ModelTree::RemoveUnref(int iBeginID, int iEndID, int iOrder)
{
	int id =  iEndID;
	while (id >= iBeginID)
	{
		//cout << id;
		if (accumulate(mModelTree[id].reference_count.begin(),mModelTree[id].reference_count.end(),0) == 0)
		{
			//cout << " Removed" << endl;
			//Decreasing reference count of arguments model tree
			// First argument is a temporary result,
			if (mModelTree[id].type1 == eTempResult)
			{					
				//Decreasing reference count of argument 1 in model tree
				int arg = mModelTree[id].id1;
				mModelTree[arg].reference_count[iOrder]--;
			}
			// Second argument has id >=0 (temporary result),
			if (mModelTree[id].id2 >= 0)
			{
				//Decreasing reference count of argument 2 in model tree
				int arg = mModelTree[id].id2;
				mModelTree[arg].reference_count[iOrder]--;
			}
			//Updating equals ids in mDerivativeIndex
			for (int d=0; d<mDerivativeIndex[iOrder-1].size();d++)
			{
				if (mDerivativeIndex[iOrder-1][d].token_id>id)
					mDerivativeIndex[iOrder-1][d].token_id--;
			}
			//cout << "ModelTree size : " << mModelTree.size() << endl;
			//Updatting upper token ids in model tree and map
			mIndexOfTokens.erase(Key((MToken) mModelTree[id]));
			for (int id2 = id+1; id2 <= iEndID; id2++)
			{
				//cout << " - " << mIndexOfTokens[Key((MToken) mModelTree[id2])];
				mIndexOfTokens.erase(Key((MToken) mModelTree[id2]));
			}
			//cout << endl;
			for (int id2 = id+1; id2 <= iEndID; id2++)
			{
				// Updating derivative ids
				map<int, int, less<int> >::iterator it;
				for (it = mModelTree[id2].p1.begin(); it != mModelTree[id2].p1.end(); it++)
				{
					int p1 = (*it).second;
					int var = (*it).first;
					//cout << "===========" << mModelTree[p1].d1[var] << "/";
					mModelTree[p1].d1[var] = id2-1;
					//cout <<  mModelTree[p1].d1[var] << endl;						
				}
				// Updating ModelTree map
				if (mModelTree[id2].type1 == eTempResult)
				{	
					if (mModelTree[id2].id1>id)
					{
						mModelTree[id2].id1--;
					}
				}
				if (mModelTree[id2].id2>id)
				{
					mModelTree[id2].id2--;
				}
			}
			mModelTree.erase(mModelTree.begin()+id);
			for (int id2 = id; id2 < iEndID; id2++)
			{
				mIndexOfTokens[Key((MToken) mModelTree[id2])] = id2;
				//cout << " - " << mIndexOfTokens[Key((MToken) mModelTree[id2])];
			}
			//Removing token from model tree
			//cout << "ModelTree size : " << mModelTree.size() << endl;


			iEndID--;
			id--;
		}
		else
		{
			id--;
			//cout << endl;
		}
	}

}
//------------------------------------------------------------------------------
*/
/*
void ModelTree::DecrementUnref(int iBeginID, int iEndID, int iOrder)
{
	int id =  iEndID;
	while (id >= iBeginID)
	{
		//cout << id;
		if (accumulate(mModelTree[id].reference_count.begin(),mModelTree[id].reference_count.end(),0) == 0)
		{
			//Decreasing reference count of arguments model tree
			// First argument is a temporary result,
			if (mModelTree[id].type1 == eTempResult)
			{					
				//Decreasing reference count of argument 1 in model tree
				int arg = mModelTree[id].id1;
				mModelTree[arg].reference_count[iOrder]--;
			}
			// Second argument has id >=0 (temporary result),
			if (mModelTree[id].id2 >= 0)
			{
				//Decreasing reference count of argument 2 in model tree
				int arg = mModelTree[id].id2;
				mModelTree[arg].reference_count[iOrder]--;
			}
			id--;
		}
		else
		{
			id--;
			//cout << endl;
		}
	}

}
*/
//------------------------------------------------------------------------------
inline string ModelTree::getArgument(NodeID id, Type type, EquationType iEquationType)
{
	
	stringstream 	argument;
	

	if (type == eParameter)
	{
		argument << param_name << lpar << (int)id+offset << rpar;
	}
	else if (type == eNumericalConstant)
	{		
		argument << NumericalConstants::get((int) id);
	}	
	else if (type == eEndogenous || type == eExogenous || type == eExogenousDet)
	  if (iEquationType == eStaticEquations || iEquationType == eStaticDerivatives)
	  {		
		int idx = VariableTable::getSymbolID((int) id)+offset;
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
			argument <<  "exedet_" << lpar << idx << rpar;
		}
	  }
	  else
	  {
		if (type == eEndogenous)
		{
			int idx = VariableTable::getPrintIndex((int) id)+offset;
			argument <<  "y" << lpar << idx << rpar;
		}
		else if (type == eExogenous)
		{	
			int idx = VariableTable::getSymbolID((int) id)+offset;
			int lag = VariableTable::getLag((int) id);
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
			int idx = VariableTable::getSymbolID((int) id)+offset;
			int lag = VariableTable::getLag((int) id);
			if (offset == 1)
			{
				if (lag != 0)
				{
					argument <<  "exdet_" << lpar << "it_ + " << lag 
						<< ", " << idx << rpar;
				}
				else
				{
					argument <<  "exdet_" << lpar << "it_, " << idx << rpar;
				}
			}   
			else
			{
				if (lag != 0)
				{
					argument <<  "exdet_" << lpar << "it_ + " << lag 
						<< "+" << idx <<  "*nb_row_xd" << rpar;
				}
				else
				{
					argument <<  "exdet_" << lpar << "it_+" << idx << "*nb_row_xd" <<  rpar;
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
		param_name = "M_.params";
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
		for (int lag = -ModelParameters::max_lag; lag <= ModelParameters::max_lead; lag++)
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
	output << "M_.exo_name_orig_ord = [1:" << ModelParameters::exo_nbr << "];\n";
	output << "M_.maximum_lag = " << ModelParameters::max_lag << ";\n";
	output << "M_.maximum_lead = " << ModelParameters::max_lead<< ";\n";
	if (ModelParameters::endo_nbr)
		output << "oo_.steady_state = zeros(" << ModelParameters::endo_nbr << ", 1);\n";
	if (ModelParameters::exo_nbr)
		output << "oo_.exo_steady_state = zeros(" << ModelParameters::exo_nbr << ", 1);\n";
	if (ModelParameters::parameter_nbr)
		output << "M_.param = zeros(" << ModelParameters::parameter_nbr << ", 1);\n";
	if (ModelParameters::exo_det_nbr)
		output << "oo_exdet_ = zeros(" << ModelParameters::exo_det_nbr << ", 1);\n";
	if (ModelParameters::exo_det_nbr)
		output << "oo_exedet_ = zeros(" << ModelParameters::exo_det_nbr << ", 1);\n";
}
//------------------------------------------------------------------------------
string ModelTree::get()
{
	return output.str();
}
//------------------------------------------------------------------------------
#ifdef TEST_MODELTREE
int main(void)
{
	SymbolTable 		st; 
	VariableTable 		vt;
	NumericalConstants 	nc;
	ModelTree			model;
	vector<int>			t(20);
	
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
	t[14] =	model.AddToken(Argument(t[13], eTempResult), 
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
	//	cout << "error---------------------\n";
	//	exit(-1);
	//}
	//cout << model.getDynamicDeterministicModel() << endl;
	//cout << model.getDynamicStochasticModel() << endl;	
	//VariableTable::Sort();
}
#endif
//------------------------------------------------------------------------------
