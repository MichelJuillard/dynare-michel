/*! \file
  \version 1.0
  \date 04/09/2004
  \par This file implements the DataTree class methodes.
*/
//------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <map>
using namespace std;
#include <time.h>

//------------------------------------------------------------------------------
#include "DynareBison.hh"
#include "VariableTable.hh"
#include "NumericalConstants.hh"
#include "DataTree.hh"

DataTree::DataTree(SymbolTable &symbol_table_arg) :
  symbol_table(symbol_table_arg),
  variable_table(symbol_table_arg),
  NoOpCode(-1), NullID(NULL)
{
  offset = 1;
  current_order = 0;

  Zero = new MetaToken(reinterpret_cast <NodeID> (0), eNumericalConstant, NULL, -1);
  Zero->op_name = "";
  Zero->reference_count.resize(current_order+1,2);
  Zero->idx = 0;
  mModelTree.push_back(Zero);
  mIndexOfTokens[*Zero]=Zero;

  One = new MetaToken(reinterpret_cast <NodeID> (1), eNumericalConstant, NULL, -1);  One->op_name = "";
  One->op_name = "";
  One->reference_count.resize(current_order+1,1);
  One->idx = 1;
  mModelTree.push_back(One);
  mIndexOfTokens[*One]=One;

  MinusOne = new MetaToken(One, eTempResult, NULL, token::UMINUS);
  MinusOne->op_name = OperatorTable::str(token::UMINUS);
  MinusOne->reference_count.resize(current_order+1,1);
  MinusOne->idx = 2;
  mModelTree.push_back(MinusOne);
  mIndexOfTokens[*MinusOne]=MinusOne;

  // Pushing "0=0" into mModelTree
  ZeroEqZero = new MetaToken(Zero, eTempResult, Zero, token::EQUAL);
  ZeroEqZero->op_name = OperatorTable::str(token::EQUAL);
  ZeroEqZero->reference_count.resize(current_order+1,1);
  ZeroEqZero->idx = 3;
  mModelTree.push_back(ZeroEqZero);
  mIndexOfTokens[*ZeroEqZero]=ZeroEqZero;

  // Initialise global node counter
  nodeCounter = 4;

  BeginModel = mModelTree.end();
  BeginModel--;
}

//------------------------------------------------------------------------------
DataTree::~DataTree()
{
  for (TreeIterator it = mModelTree.begin(); it != mModelTree.end(); it++)
    {
      if (*it != NullID)  delete *it;
    }
}
