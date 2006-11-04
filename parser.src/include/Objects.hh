#ifndef _OBJECTS_HH
#define _OBJECTS_HH
//------------------------------------------------------------------------------
/*! \file
  \version 1.0
  \date 04/13/2004
  \par This file defines the Objects class.
*/
//------------------------------------------------------------------------------
using namespace std;
#include "SymbolTable.hh"
#include "ModelTypes.hh"
//------------------------------------------------------------------------------
namespace dynare
{
  /*!
    \class  Objects
    \brief  This class defines data associated to parsed tokens.
  */
  class Objects
  {
  public :
    /*! Parsed string name : can be literal name or operator */
    string      symbol;
    /*! ID of object : can be ID of symbol, variable, or token in model tree */
    NodeID      ID;
    /*! Type of object */
    Type      type;
    /*! In case of operator object, this is set to its code */
    int       opcode;
  public :
    /*! Constructor with default values */
    Objects()
    {
      ID = NULL;
      symbol = "";
      type = eUNDEF;
    }
    /*! Constructor of object with known symbol name */
    Objects(string name)
    {
      opcode = NAME;
      symbol = name;
    }
    /*! Constructor of object with known symbol name */
    Objects(const char* name, NodeID id = NULL, Type t = eUNDEF)
    {
      opcode = NAME;
      symbol = name;
      type = t;
      ID = id;
    }
    /*! Constructor of object with known symbol name, ID and type */
    Objects(string name,NodeID id, Type t)
    {
      symbol = name;
      type = t;
      ID = id;
    }
    /*! Constructor of object with known operator code */
    Objects(int op)
    {
      opcode = op;
    }
    /*! conversion to string */
    operator string() const
    {
      return symbol;
    }
    /*! Copy constructor */
    Objects(const Objects& obj) : symbol(obj.symbol)
    {
      // Empty
    }
    /*! Destructor */
    ~Objects()
    {
      // Empty
    };
  };

}
#endif
