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
//#include <ext/hash_map>
using namespace std;
//using __gnu_cxx::hash;
//using __gnu_cxx::hash_map;
#include <time.h>

//------------------------------------------------------------------------------
#include "DynareBison.hh"
#include "VariableTable.hh"
#include "NumericalConstants.hh"
#include "DataTree.hh"
//------------------------------------------------------------------------------
const int DataTree::NoOpCode = -1;
const NodeID DataTree::NullID = NULL;
const NodeID DataTree::Zero = new MetaToken(reinterpret_cast <NodeID> (0), eNumericalConstant, NULL, -1);
const NodeID DataTree::One = new MetaToken(reinterpret_cast <NodeID> (1), eNumericalConstant, NULL, -1);
const NodeID DataTree::MinusOne = new MetaToken(One, eTempResult, NULL, UMINUS);
const NodeID DataTree::ZeroEqZero = new MetaToken(Zero, eTempResult, Zero, EQUAL);
int      DataTree::offset = 1;
//------------------------------------------------------------------------------
DataTree::DataTree()
{

  current_order = 0;
  //Here "0" and "1" have been added to NumericalConstants class
  SymbolTable::AddSymbolDeclar("0.0",eNumericalConstant, "");
  SymbolTable::AddSymbolDeclar("1.0",eNumericalConstant, "");

  Zero->op_name = "";
  Zero->reference_count.resize(current_order+1,2);
  Zero->idx = 0;
  mModelTree.push_back(Zero);
  mIndexOfTokens[Zero->Key()]=Zero;

  One->op_name = "";
  One->reference_count.resize(current_order+1,1);
  One->idx = 1;
  mModelTree.push_back(One);
  mIndexOfTokens[One->Key()]=One;

  MinusOne->op_name = operator_table.str(UMINUS);
  MinusOne->reference_count.resize(current_order+1,1);
  MinusOne->idx = 2;
  mModelTree.push_back(MinusOne);
  mIndexOfTokens[MinusOne->Key()]=MinusOne;

  // Pushing "0=0" into mModelTree
  ZeroEqZero->op_name = operator_table.str(EQUAL);
  ZeroEqZero->reference_count.resize(current_order+1,1);
  ZeroEqZero->idx = 3;
  mModelTree.push_back(ZeroEqZero);
  mIndexOfTokens[ZeroEqZero->Key()]=ZeroEqZero;

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
