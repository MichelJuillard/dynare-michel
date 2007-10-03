#ifndef _SYMBOLTABLETYPES_HH
#define _SYMBOLTABLETYPES_HH

enum UnaryOpcode
  {
    oUminus,
    oExp,
    oLog,
    oLog10,
    oCos,
    oSin,
    oTan,
    oAcos,
    oAsin,
    oAtan,
    oCosh,
    oSinh,
    oTanh,
    oAcosh,
    oAsinh,
    oAtanh,
    oSqrt
  };

enum BinaryOpcode
  {
    oPlus,
    oMinus,
    oTimes,
    oDivide,
    oPower,
    oEqual
  };


//! Symbol type enum
enum Type
  {
    eEndogenous = 0,               //!< Endogenous
    eExogenous = 1,                //!< Exogenous
    eExogenousDet = 2,             //!< Exogenous deterministic (new)
    eRecursiveVariable = 3,        //!< Recursive variable (reserved for future use)
    eParameter = 4,                //!< Parameter
    eModelLocalVariable = 10,      //!< Local variable whose scope is model (pound expression)
    eModFileLocalVariable = 11     //!< Local variable whose scope is mod file (model excluded)
  };

struct Symbol
{
  //! Symbol type
  Type type;
  //! Symbol ID : for each type
  int id;

  Symbol() : id(-1)
  {
  }
};

#endif
