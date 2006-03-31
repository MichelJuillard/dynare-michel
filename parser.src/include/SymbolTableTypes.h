#ifndef SYMBOLTABLETYPES_H
#define SYMBOLTABLETYPES_H
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
	eEndogenous = 0,		//!< Endogenousous 
	eExogenous = 1,			//!< Exogenousous 
	eExogenousDet = 2, 		//!< Exogenousous deterministic (new) 
	eRecursiveVariable = 3,	//!< Recursive variable (reserved for future use) 
	eParameter = 4,			//!< Parameter 
	eLocalParameter = 41,           //!< Parameter  local to a model
	eLoopIndex = 5,			//!< Loop index 
	eTempResult = 6,		//!< Temporary result, used only in Expression class 
	eNumericalConstant = 7,	//!< Numerical constant,  used only in Expression class  
	eUnkownFunction = 8,	//!< Unkown functions, used only in Expression class  
	eUNDEF = 9				//!< Undefinded 
};
/*! Symbol reference flag enum */
enum Reference
{
	eNotReferenced,		//!< Not yet referenced in model 
	eReferenced,		//!< Already referenced in model 
	eOutOfScoop,		//!< Out of scoop (for loop index) 
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
