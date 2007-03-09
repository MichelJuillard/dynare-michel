#ifndef _SYMBOLTABLETYPES_HH
#define _SYMBOLTABLETYPES_HH

//! Symbol type enum
enum Type
  {
    eEndogenous = 0,               //!< Endogenous
    eExogenous = 1,                //!< Exogenous
    eExogenousDet = 2,             //!< Exogenous deterministic (new)
    eRecursiveVariable = 3,        //!< Recursive variable (reserved for future use)
    eParameter = 4,                //!< Parameter
    eLocalParameter = 10,          //!< Parameter  local to a model
  };

//! Symbol reference flag enum
enum Reference
  {
    eNotReferenced,                //!< Not yet referenced in model
    eReferenced,                   //!< Already referenced in model
  };

struct Symbol
{
  //! Symbol type
  Type type;
  //! Symbol ID : for each type
  int id;
  //! Symbol reference flag
  Reference referenced;

  Symbol() : id(-1), referenced(eNotReferenced)
  {
  }
};

#endif
