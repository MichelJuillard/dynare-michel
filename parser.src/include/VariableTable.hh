#ifndef _VARIABLETABLE_HH
#define _VARIABLETABLE_HH
//------------------------------------------------------------------------------
/** \file 
 * \version 1.0
 * \date 12/16/2003
 * \par This file defines the VariableTable class .
 */
//------------------------------------------------------------------------------
#include <map>
#include <string>
#include <vector>
//------------------------------------------------------------------------------
#include "SymbolTable.hh"
//------------------------------------------------------------------------------
/*! 
 \class Variable
 \brief	Variable struct 
*/
struct Variable {
  /*! Variable type */
  Type 	type;
  /*! Symbol ID */
  int  	symbol_id;
  /*! Variable ID */
  int 	variable_id;
};
/*! Variable key type to acced variable table elements */
typedef std::pair<std::string, int> varKey;
//------------------------------------------------------------------------------
/*!
 \class VariableTable 
 \brief  This class is used to store variables as they appear 
 in the model (with their lead or lag)
*/
class VariableTable
{ 
	private :
		/*!  Variable table data */
		//static std::map<varKey,Variable>	mVariableTable;
		static std::map<varKey,int>	mVariableTable;
		/*! Index (IDs) of variables in variable table */
		static std::vector<varKey>			mVariableIndex;
		/*! Variable IDs of sorted variable table */
		static std::vector<int>			mSortedVariableID;
		/*! Output index for variable table */
		static std::vector<int>			mPrintFormatIndex;			
	public :
		/*! */
		VariableTable();
		/*! */
		~VariableTable();
		/*! Find type and ID in SymbolTable 
		 - Increment variable_id; 
		 - Make variable 
		 - Push variable on variabletable 
		*/	
		static int 	AddVariable(std::string iName, int iLag);	
		/*! Pointer to error function of parser class */                               
		static void (* error) (const char* m);
		/*! Decremente a symbol id of a variable */
		static void decSymbolID(std::string iName, int id, int iLag, Type iType);
		/*! Return VariableTable[name,lag].variable_id  */
		inline static int 	getID(std::string iName, int iLag);
		/*! Return lag of variable */
		inline static int 	getLag(int iID);
		/*! Return symbol ID of variable */
		inline static int  getSymbolID(int ivarID);
		/*! Gets varibale type */
		inline static Type	getType(int ivarID);
		/*! Gets nomber of variables in mVariableTable */
		inline static int  size(void);
		/*! Gets variable ID of sorted variable table */
		inline static int  getSortID(int);
		/*! Return variable index to print in format : y(index) or oo_.y_simul(index) ... */
		inline static int 	getPrintIndex(int iVarID);
		/*! Sorts variable table */
		static void	Sort(void);
};
inline int  VariableTable::getSortID(int iVarID)
{
	return mSortedVariableID[iVarID];
}
//------------------------------------------------------------------------------
inline int VariableTable::getPrintIndex(int iVarID)
{
	return mPrintFormatIndex[iVarID];
}
//------------------------------------------------------------------------------
inline int VariableTable::getID(std::string iName, int iLag)
{
	if (mVariableTable.find(make_pair(iName, iLag)) == mVariableTable.end())
	{
		return -1;
	}
	else
	{
		return mVariableTable[make_pair(iName, iLag)];
	}
}
//------------------------------------------------------------------------------
inline Type VariableTable::getType(int ivarID)
{
	varKey key = mVariableIndex[ivarID];
	//return mVariableTable[key].type;
	return SymbolTable::getType(key.first);
}
//------------------------------------------------------------------------------
inline int VariableTable::getSymbolID(int ivarID)
{
	varKey key = mVariableIndex[ivarID];
	//return mVariableTable[key].symbol_id;
	return SymbolTable::getID(key.first);
}
//------------------------------------------------------------------------------
inline int VariableTable::getLag(int iID)
{
	return mVariableIndex[iID].second;
}
//------------------------------------------------------------------------------
inline int VariableTable::size(void)
{
	return mVariableTable.size();
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
#endif
