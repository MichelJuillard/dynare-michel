#ifndef _SYMBOLTABLETYPES_HH
#define _SYMBOLTABLETYPES_HH
//------------------------------------------------------------------------------
/*! \file
  \version 1.0
  \date 04/26/2004
  \par This file defines types related to SymbolTable.
*/
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/*! Symbol type enum */
enum Type
  {
    eEndogenous = 0,               //!< Endogenous
    eExogenous = 1,                //!< Exogenous
    eExogenousDet = 2,             //!< Exogenous deterministic (new)
    eRecursiveVariable = 3,        //!< Recursive variable (reserved for future use)
    eParameter = 4,                //!< Parameter
    eLocalParameter = 10,          //!< Parameter  local to a model
    eLoopIndex = 5,                //!< Loop index
    eTempResult = 6,               //!< Temporary result, used only in Expression class
    eNumericalConstant = 7,        //!< Numerical constant,  used only in Expression class
    eUNDEF = 9                     //!< Undefined
  };
/*! Symbol reference flag enum */
enum Reference
  {
    eNotReferenced,                //!< Not yet referenced in model
    eReferenced,                   //!< Already referenced in model
  };
/*!
  \class Symbol
  \brief Symbol structure
*/
struct Symbol
{
  /*! Symbol type */
  Type type;
  /*! Symbol ID : for each type */
  int id;
  /*! Symbol reference flag */
  Reference referenced;
  /*! Lags of symbol if it is a variable */
  std::vector<int> lags;
  Symbol()
  {
    type = eUNDEF;
    id = -1;
    referenced = eNotReferenced;
  }
} ;
//------------------------------------------------------------------------------
#endif
