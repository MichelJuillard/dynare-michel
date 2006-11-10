#ifndef _DATATREE_HH
#define _DATATREE_HH
//------------------------------------------------------------------------------
/*! \file
  \version 1.0
  \date 04/09/2004
  \par This file defines the DataTree class.
*/
//------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <list>
#include <stack>
#include <sstream>
#include <fstream>
#include <map>
#include <stdio.h>
//------------------------------------------------------------------------------
#include "SymbolTable.hh"
#include "OperatorTable.hh"
#include "NumericalConstants.hh"
#include "ModelTypes.hh"
#include "VariableTable.hh"
//------------------------------------------------------------------------------

typedef std::map<MToken, NodeID, MTokenLess> TreeMap;
typedef std::list<NodeID> TreeList;
typedef TreeList::iterator TreeIterator;
/*!
  \class  DataTree
  \brief  Provides data storage for model tree
*/
class DataTree
{
protected :
  /*! A list of structures "token" */
  TreeList    mModelTree;
  /*! matches key with entry id in list of token */
  TreeMap         mIndexOfTokens;
  /*! A counter for filling MetatToken's idx field */
  int nodeCounter;
  /*! Operator table : names and precedence of oparator */
  OperatorTable     operator_table;
  /*! ID of first token in model tree (first tokens are "0", "1" and "0=" */
  TreeIterator      BeginModel;
  /*!
    Used with field reference_count to know if
    dealing with derivative 0, 1, 2 etc
  */
  int         current_order;
  /*! Pushs token into model tree */
  inline NodeID     PushToken(NodeID iArg1,int iOpCode, NodeID iArg2 = NULL, Type iType1 = eTempResult);
  inline NodeID   getIDOfToken(const MToken &iToken);
public :
  /*! Flag for empty operator (final toekn) */
  static const int        NoOpCode;
  static const NodeID     NullID;
  static const NodeID     Zero;
  static const NodeID     One;
  static const NodeID     MinusOne;
  static const NodeID     ZeroEqZero;
  /*! Type of output 0 for C and 1 for Matlab (default) , also used as matrix index offset*/
  static int    offset;
  /*! Pointer to error function of parser class */
  void (* error) (const char* m);
  /*! Increment reference count of given token */
  inline void     IncrementReferenceCount(NodeID token);
  /*! Adds terminal token to model tree */
  inline NodeID     AddTerminal(NodeID iArg, Type type);
  /*! Adds terminal  token (a variable) to model tree */
  inline NodeID     AddTerminal(std::string iArgName, int iLag = 0);
  /*! Adds "arg1+arg2" to model tree */
  inline NodeID     AddPlus(NodeID iArg1, NodeID iArg2);
  /*! Adds "arg1-arg2" to model tree */
  inline NodeID     AddMinus(NodeID iArg1, NodeID iArg2);
  /*! Adds "-arg" to model tree */
  inline NodeID     AddUMinus(NodeID iArg1);
  /*! Adds "arg1*arg2" to model tree */
  inline NodeID     AddTimes(NodeID iArg1, NodeID iArg2);
  /*! Adds "arg1/arg2" to model tree */
  inline NodeID     AddDivide(NodeID iArg1, NodeID iArg2);
  /*! Adds "arg1^arg2" to model tree */
  inline NodeID     AddPower(NodeID iArg1, NodeID iArg2);
  /*! Adds "exp(arg)" to model tree */
  inline NodeID     AddExp(NodeID iArg1);
  /*! Adds "log(arg)" to model tree */
  inline NodeID     AddLog(NodeID iArg1);
  /*! Adds "log10(arg)" to model tree */
  inline NodeID     AddLog10(NodeID iArg1);
  /*! Adds "cos(arg)" to model tree */
  inline NodeID     AddCos(NodeID iArg1);
  /*! Adds "sin(arg)" to model tree */
  inline NodeID     AddSin(NodeID iArg1);
  /*! Adds "tan(arg)" to model tree */
  inline NodeID     AddTan(NodeID iArg1);
  /*! Adds "acos(arg)" to model tree */
  inline NodeID     AddACos(NodeID iArg1);
  /*! Adds "asin(arg)" to model tree */
  inline NodeID     AddASin(NodeID iArg1);
  /*! Adds "atan(arg)" to model tree */
  inline NodeID     AddATan(NodeID iArg1);
  /*! Adds "cosh(arg)" to model tree */
  inline NodeID     AddCosH(NodeID iArg1);
  /*! Adds "sinh(arg)" to model tree */
  inline NodeID     AddSinH(NodeID iArg1);
  /*! Adds "tanh(arg)" to model tree */
  inline NodeID     AddTanH(NodeID iArg1);
  /*! Adds "acosh(arg)" to model tree */
  inline NodeID     AddACosH(NodeID iArg1);
  /*! Adds "asinh(arg)" to model tree */
  inline NodeID     AddASinH(NodeID iArg1);
  /*! Adds "atanh(args)" to model tree */
  inline NodeID     AddATanH(NodeID iArg1);
  /*! Adds "sqrt(arg)" to model tree */
  inline NodeID     AddSqRt(NodeID iArg1);
  /*! Adds "arg1=arg2" to model tree */
  inline NodeID     AddEqual(NodeID iArg1, NodeID iArg2);
  /*! Adds "arg1=arg2" as assignment to model tree */
  inline NodeID     AddAssign(NodeID iArg1, NodeID iArg2);
public :
  /*! Constructor */
  DataTree();
  /*! Destructor */
  ~DataTree();
};
//------------------------------------------------------------------------------
inline NodeID DataTree::PushToken(NodeID iArg1,int iOpCode, NodeID iArg2, Type iType1)
{
  MetaToken*  lToken = new MetaToken(iArg1, iType1, iArg2, iOpCode);

  if (iOpCode != NoOpCode)
    {
      lToken->op_name = operator_table.str(iOpCode);
    }
  else
    {
      lToken->op_name = "";
    }
  lToken->reference_count.resize(current_order+1,0);
  lToken->idx = nodeCounter++;
  mModelTree.push_back(lToken);
  //Updating reference counters and time costs
  if (iType1 == eTempResult )
    {
      lToken->cost = operator_table.cost(iOpCode,offset)+iArg1->cost;
      IncrementReferenceCount(iArg1);
    }
  if(iArg2 != NullID)
    {
      lToken->cost += iArg2->cost;
      IncrementReferenceCount(iArg2);
    }
  mIndexOfTokens[*lToken] = lToken;

  /*
    std::cout << "ID = " << ID << " / " << mIndexOfTokens.size()-1<< " - " << getIDOfToken(lToken2) << " " <<
    lToken.op_code <<  " " << lToken.id1   << " " <<
    lToken.id2 << " " << lToken.type1<<  "\n";
  */
  /*
    std::cout << "\nmIndexOfTokens------------\n";
    for (iter=mIndexOfTokens.begin();iter != mIndexOfTokens.end(); iter++)
    std::cout << (*iter).second << " -> " << (*iter).first << std::endl;
  */
  return lToken;
}

//------------------------------------------------------------------------------
inline void DataTree::IncrementReferenceCount(NodeID token)
{
  int s = token->reference_count.size();
  for (int i = s; i <= current_order; i++)
    {
      int rc = token->reference_count[i-1];
      token->reference_count.push_back(rc);
    }
  token->reference_count[current_order]++;

}

//------------------------------------------------------------------------------
inline NodeID DataTree::getIDOfToken(const MToken &iToken)
{
  TreeMap::iterator iter = mIndexOfTokens.find(iToken);
  if (iter != mIndexOfTokens.end())
    return iter->second;
  else
    return NullID;
}

//------------------------------------------------------------------------------
inline NodeID DataTree::AddTerminal(NodeID iArg1, Type iType1)
{

  if (iType1 == eTempResult)
    {
      return iArg1;
    }
  else
    {
      MToken  lToken(iArg1, iType1, NullID, NoOpCode);
      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1, NoOpCode, NullID, iType1);
    }

}

inline NodeID DataTree::AddTerminal(std::string iArgName, int iLag)
{
  int   id;
  Type  type;

  if (!SymbolTable::Exist(iArgName))
    {
      std::string msg = "unknown symbol: " + iArgName;
      (* error) (msg.c_str());
      exit(-1);
    }
  type = SymbolTable::getType(iArgName);
  if (type != eUNDEF)
    {
      SymbolTable::SetReferenced(iArgName);
      if (type == eEndogenous ||
          type == eExogenousDet  ||
          type == eExogenous  ||
          type == eRecursiveVariable)
        {
          id = VariableTable::getID(iArgName,iLag);
          if (id == -1)
            {
              std::string msg = "unknown variable " + iArgName;
              (* error) (msg.c_str());
              exit(-1);
            }
        }
      else
        id = SymbolTable::getID(iArgName);

    }
  else
    {
      std::string msg = "unknown parameter: " + iArgName;
      (* error) (msg.c_str());
      exit(-1);
    }
  return AddTerminal((NodeID) id, type);
}

inline NodeID DataTree::AddPlus(NodeID iArg1, NodeID iArg2)
{

  if (iArg1 != Zero && iArg2 != Zero)
    {
      MToken  lToken;

      if (iArg1 <= iArg2)
        lToken = MToken(iArg1, eTempResult, iArg2, PLUS);
      else
        lToken = MToken(iArg2, eTempResult, iArg1, PLUS);
      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }

      // To treat commutativity of "+"
      // Token id (iArg1 and iArg2) are sorted
      if (iArg1 <= iArg2)
        return PushToken(iArg1,PLUS,iArg2);
      else
        return PushToken(iArg2,PLUS,iArg1);
    }
  else if (iArg1 != Zero)
    {
      return iArg1;
    }
  else if (iArg2 != Zero)
    {
      return iArg2;
    }
  else
    {
      return Zero;
    }
}

inline NodeID DataTree::AddMinus(NodeID iArg1, NodeID iArg2)
{

  if (iArg1 != Zero && iArg2 != Zero)
    {
      MToken  lToken(iArg1, eTempResult, iArg2, MINUS);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1,MINUS, iArg2);
    }
  else if (iArg1 != Zero)
    {
      return iArg1;
    }
  else if (iArg2 != Zero)
    {
      return AddUMinus(iArg2);
    }
  else
    {
      return Zero;
    }
}

inline NodeID DataTree::AddUMinus(NodeID iArg1)
{
  if (iArg1 != Zero )
    {
      MToken  lToken(iArg1, eTempResult, NullID, UMINUS);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      if (iArg1->type1 == eTempResult &&
          iArg1->op_code == UMINUS)
        {
          IncrementReferenceCount(iArg1->id1);
          return iArg1->id1;

        }
      return PushToken(iArg1, UMINUS);
    }
  else
    {
      return Zero;
    }

}

inline NodeID DataTree::AddTimes(NodeID iArg1, NodeID iArg2)
{

  if (iArg1 != Zero && iArg1 != One && iArg2 != Zero && iArg2 != One)
    {
      MToken  lToken;

      if (iArg1 <= iArg2)
        lToken = MToken(iArg1, eTempResult, iArg2, TIMES);
      else
        lToken = MToken(iArg2, eTempResult, iArg1, TIMES);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      // Testing if arg[1 or 2] is "-1", in this case add "-arg" instead of "-1*iarg"
      // Here 2 is ID of "-1"
      if (iArg1 == MinusOne)
        return AddUMinus(iArg2);
      else if (iArg2 == MinusOne)
        return AddUMinus(iArg1);
      else
        {
          // To treat commutativity of "*"
          // Token id (iArg1 and iArg2) are sorted
          if (iArg1 <= iArg2)
            return PushToken(iArg1, TIMES, iArg2);
          else
            return PushToken(iArg2, TIMES, iArg1);
        }
    }
  else if (iArg1 != Zero && iArg1 != One && iArg2 == One)
    {
      return iArg1;
    }
  else if (iArg2 != Zero && iArg2 != One && iArg1 == One)
    {
      return iArg2;
    }
  else if (iArg2 == One && iArg1 == One)
    {
      return One;
    }
  else
    {
      return Zero;
    }
}

inline NodeID DataTree::AddDivide(NodeID iArg1, NodeID iArg2)
{

  if (iArg1 != Zero && iArg2 != Zero && iArg2 != One)
    {
      MToken  lToken(iArg1, eTempResult, iArg2, DIVIDE);
      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1,DIVIDE,iArg2);
    }
  else if (iArg2 == One)
    {
      return iArg1;
    }
  else if (iArg1 == Zero && iArg2 != Zero)
    {
      return Zero;
    }
  else
    {
      std::cout << "DIVIDE 0/0 non available\n";
      exit(-1);
    }
}

inline NodeID DataTree::AddPower(NodeID iArg1, NodeID iArg2)
{
  if (iArg1 != Zero && iArg2 != Zero && iArg2 != One)
    {
      MToken  lToken(iArg1, eTempResult, iArg2, POWER);
      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1,POWER,iArg2);
    }
  else if (iArg2 == One)
    {
      return iArg1;
    }
  else if (iArg2 == Zero)
    {
      return One;
    }
  else
    {
      return Zero;
    }

}

inline NodeID DataTree::AddExp(NodeID iArg1)
{

  if (iArg1 != Zero)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }

      return PushToken(iArg1, EXP);
    }
  else
    {
      return One;
    }
}

inline NodeID DataTree::AddLog(NodeID iArg1)
{
  if (iArg1 != Zero && iArg1 != One)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1, LOG);
    }
  else if (iArg1 == One)
    {
      return Zero;
    }
  else
    {
      std::cout << "log(0) isn't available\n";
      exit(-1);
    }
}

inline NodeID DataTree::AddLog10(NodeID iArg1)
{
  if (iArg1 != Zero && iArg1 != One)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }

      return PushToken(iArg1, LOG);
    }
  else if (iArg1 == One)
    {
      return Zero;
    }
  else
    {
      std::cout << "log10(0) isn't available\n";
      exit(-1);
    }
}

inline NodeID DataTree::AddCos(NodeID iArg1)
{
  if (iArg1 != Zero)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1, COS);
    }
  else
    {
      return One;
    }
}

inline NodeID DataTree::AddSin(NodeID iArg1)
{
  if (iArg1 != Zero)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1, SIN);
    }
  else
    {
      return Zero;
    }
}

inline NodeID DataTree::AddTan(NodeID iArg1)
{
  if (iArg1 != Zero)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1, TAN);
    }
  else
    {
      return Zero;
    }
}

inline NodeID DataTree::AddACos(NodeID iArg1)
{
  MToken  lToken(iArg1, eTempResult, NullID, ACOS);

  NodeID ID = getIDOfToken(lToken);
  if (ID != NullID)
    {
      return ID;
    }
  return PushToken(iArg1, COS);
}

inline NodeID DataTree::AddASin(NodeID iArg1)
{
  if (iArg1 != Zero)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1, SIN);
    }
  else
    {
      return Zero;
    }
}

inline NodeID DataTree::AddATan(NodeID iArg1)
{
  if (iArg1 != Zero)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1, ATAN);
    }
  else
    {
      return Zero;
    }
}

inline NodeID DataTree::AddCosH(NodeID iArg1)
{
  if (iArg1 != Zero)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1, COSH);
    }
  else
    {
      return Zero;
    }
}

inline NodeID DataTree::AddSinH(NodeID iArg1)
{
  if (iArg1 != Zero)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1, SINH);
    }
  else
    {
      return Zero;
    }
}

inline NodeID DataTree::AddTanH(NodeID iArg1)
{
  if (iArg1 != Zero)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1, TANH);
    }
  else
    {
      return Zero;
    }
}

inline NodeID DataTree::AddACosH(NodeID iArg1)
{
  MToken  lToken(iArg1, eTempResult, NullID, ACOSH);

  NodeID ID = getIDOfToken(lToken);
  if (ID != NullID)
    {
      return ID;
    }
  return PushToken(iArg1, ACOSH);
}

inline NodeID DataTree::AddASinH(NodeID iArg1)
{
  if (iArg1 != Zero)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1, ASINH);
    }
  else
    {
      return Zero;
    }
}

inline NodeID DataTree::AddATanH(NodeID iArg1)
{
  if (iArg1 != Zero)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1, ATANH);
    }
  else
    {
      return Zero;
    }
}

inline NodeID DataTree::AddSqRt(NodeID iArg1)
{
  if (iArg1 != Zero)
    {
      MToken  lToken(iArg1, eTempResult, NullID, EXP);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        {
          return ID;
        }
      return PushToken(iArg1, SQRT);
    }
  else
    {
      return Zero;
    }
}

inline NodeID DataTree::AddEqual(NodeID iArg1, NodeID iArg2=Zero)
{
  if (iArg1 == Zero && iArg2 == Zero)
    return ZeroEqZero;
  else
    {
      MToken  lToken(iArg1, eTempResult, iArg2, EQUAL);

      NodeID ID = getIDOfToken(lToken);
      if (ID != NullID)
        return ID;

      return PushToken(iArg1, EQUAL, iArg2);
    }
}

inline NodeID DataTree::AddAssign(NodeID iArg1, NodeID iArg2=Zero)
{
  MToken  lToken(iArg1, eTempResult, iArg2, ASSIGN);

  NodeID ID = getIDOfToken(lToken);
  if (ID != NullID)
    {
      return ID;
    }
  return PushToken(iArg1,ASSIGN,iArg2);
}

//------------------------------------------------------------------------------
#endif
