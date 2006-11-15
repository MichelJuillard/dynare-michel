#ifndef _MODELTYPES_HH
#define _MODELTYPES_HH
//------------------------------------------------------------------------------
/*! \file
  \version 1.0
  \date 04/13/2004
  \par This file defines types for ModelTree class.
*/
//------------------------------------------------------------------------------
#include <string>
#include <map>

typedef class MetaToken* NodeID;

/*!
  \struct  DerivativeIndex
  \brief  Derivative index structure
*/
//------------------------------------------------------------------------------
struct DerivativeIndex
{
  /*!/ id of "equal token" in mModelTree */
  NodeID token_id;
  /*! id of equation being derived */
  int equation_id;
  /*! column number of Jacobian or Hessian matrix */
  int derivators;
  /*! Constructor */
  inline DerivativeIndex(NodeID iToken_id, int iEquationID, int iDerivator)
  {
    token_id = iToken_id;
    equation_id = iEquationID;
    derivators = iDerivator;
  }
};
/*!
  \struct  MToken
  \brief  Basic token structure, essentially for computing hash keys
*/
//------------------------------------------------------------------------------
struct MToken
{
public :
  /*! ID of first operand */
  NodeID  id1;
  /*! Type of first operand */
  Type  type1;
  /*! ID of second operand */
  NodeID  id2;
  /*! Operator code */
  int   op_code;
  /*! Constructor */
  inline MToken(NodeID iId1 = NULL, Type iType1 = eUNDEF, NodeID iId2 = NULL, int iOpCode = -1)
  {
    id1 = iId1;
    type1 = iType1;
    id2 = iId2;
    op_code = iOpCode;

  }
  /*! Copy constructor */
  inline MToken(const MToken& t)
  {
    id1 = t.id1;
    id2 = t.id2;
    type1 = t.type1;
    op_code = t.op_code;
  }
  /*! Destructor */
  ~MToken() { }
};

/*!
  \struct MTokenLess
  \brief Class which defines a lexicographic order over MToken, used in std::map
*/
struct MTokenLess
{
  bool operator()(const MToken &n1, const MToken &n2) const
  {
    if (n1.id1 != n2.id1)
      return(n1.id1 < n2.id1);
    else if (n1.id2 != n2.id2)
      return(n1.id2 < n2.id2);
    else if (n1.type1 != n2.type1)
      return(n1.type1 < n2.type1);
    else
      return(n1.op_code < n2.op_code);
  }
};

//------------------------------------------------------------------------------
/*!
  \struct  MetaToken
  \brief  Meta token structure
*/
struct MetaToken : public MToken
{
private:
  /*! Address of first order partial derivative */
  std::map<int, NodeID, std::less<int> >  d1;
public :
  /*! Index of token starting from zero */
  int             idx;
  /*!
    Number of times a given token is referenced?
    Each element of vector refer to derivative 0, 1, 2, etc
  */
  std::vector<int>        reference_count;
  /*! Commulative computation time cost */
  int               cost;
  /*!
    If the token can be a temporary result, this flag is set to 1.
    If it can not be a temporary result (operation with constants etc)
    it is set to 0.
    After writing corresponding expression for the first time it is set to -1.
  */
  int             tmp_status;
  /*! Name of unkown functions */
  std::string             func_name;
  /*! Expression of token, exp[0] for static model
    exp[1] for dynamic model */
  std::string           exp[2];
  /*! Flag : operand 1 is followed or not */
  bool            followed1;
  /*! Flag : operand 2 is followed or not */
  bool            followed2;
  /*! Operator name */
  std::string           op_name;
  /*! set to one when left node has been treated */
  int left_done;
  /*! set to one when right node has been treated */
  int right_done;
  /*! set to one if closing parenthesis after token */
  int close_parenthesis;

  /*! Ordered list of variable IDs with respect to which the derivative is potentially non-null */ 
  std::vector<int> non_null_derivatives;

  /*! Constructor */
  inline MetaToken(NodeID iId1 = NULL, Type iType1 = eUNDEF,NodeID iId2= NULL, int iOpCode = -1)
    : MToken(iId1, iType1,iId2,iOpCode)
  {
    cost = 0;
    followed1 = followed2 = false;
    tmp_status = 0;
    left_done = 0;
    right_done = 0;
    close_parenthesis = 0;

    switch(type1)
      {
      case eEndogenous:
      case eExogenous:
      case eExogenousDet:
      case eRecursiveVariable:
        // For a variable, the only non-null derivative is with respect to itself
        non_null_derivatives.push_back((long int) id1);
        break;
      case eParameter:
      case eLocalParameter:
      case eNumericalConstant:
      case eUNDEF:
        // All derivatives are null
        break;
      case eTempResult:
        if (!id2)
          non_null_derivatives = id1->non_null_derivatives;
        else
          {
            // Compute set union of id1->non_null_derivates and id2->non_null_derivatives
            // The result is ordered by ascending order
            std::vector<int>::iterator it1, it2;
            it1 = id1->non_null_derivatives.begin();
            it2 = id2->non_null_derivatives.begin();
            while(it1 != id1->non_null_derivatives.end()
                  && it2 != id2->non_null_derivatives.end())
              {
                if (*it1 == *it2)
                  {
                    non_null_derivatives.push_back(*it1++);
                    it2++;
                  }
                else if (*it1 < *it2)
                  non_null_derivatives.push_back(*it1++);
                else
                  non_null_derivatives.push_back(*it2++);
              }
            while(it1 != id1->non_null_derivatives.end())
                non_null_derivatives.push_back(*it1++);
            while(it2 != id2->non_null_derivatives.end())
              non_null_derivatives.push_back(*it2++);
          }
        break;
      }
  }
  /*! Destructor */
  ~MetaToken()
  {
  }

  /*! Sets derivative with respect to iVarID equal to iDerivative */
  inline void setDerivativeAddress(NodeID iDerivative, int iVarID)
  {
    d1[iVarID] = iDerivative;
  }

  /*! Get derivative with respect to iVarID.
    Defined in ModelTree.cc because it needs DataTree::Zero and DataTree::ZeroEqZero
    and is only used in that source file.
  */
  inline NodeID getDerivativeAddress(int iVarID);
};
//------------------------------------------------------------------------------
/*! Equation type enum */
enum EquationType
  {
    eStaticEquations = 0,          //!< Equation of static model
    eDynamicEquations = 1,         //!<  Equation of dynamic model
    eStaticDerivatives = 2,        //!< Derivative of static model
    eDynamicDerivatives = 3        //!< Derivative of dynamic model
  };
//------------------------------------------------------------------------------
#endif
