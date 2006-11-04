#ifndef _OPERATORTABLE_HH
#define _OPERATORTABLE_HH
//------------------------------------------------------------------------------
/** \file 
 * \version 1.0
 * \date 04/26/2004
 * \par This file defines the OperatorTable class .
 */
//------------------------------------------------------------------------------
#include <map>
#include <string>
#include <vector>
//------------------------------------------------------------------------------
#include "DynareBison.hh"
//------------------------------------------------------------------------------
/*!
  \class  OperatorTable
  \brief  Defines operators used in ModelTree
*/
//------------------------------------------------------------------------------
class OperatorTable
{
	private :
		/*! 
		 \class Operator
		 \brief Operator structure 
		 */
		struct Operator
		{
			Operator()
			{
				cost.resize(2);
			}
			/*! Operator name for Matlab */
			std::string			str;
			/*! Operator precedence */
			int				precedence;
			/*! Set to true if operator is a function */
			bool			isfunction;
			/*! Time computation cost of operator */
			std::vector<int>				cost;
		};
		/*! Type of operators map indexed by code */
		typedef std::map<int, Operator> operator_map;
		/*! Operator table */
		operator_map	operator_table;
	public :
		/*! Constructor */
		OperatorTable();
		/*! Destrcutor */
		~OperatorTable();
		/*! Sets oparator name*/
		inline std::string 			str(int op_code);
		/*! Sets operator precedence */
		inline int				precedence(int op_code);
		/*! Sets flag to true if operator is a function */
		inline bool			isfunction(int op_code);
		/*! Gets cost of operator */
		inline int 			cost(int op_code, int what);
};
//------------------------------------------------------------------------------
inline std::string 	OperatorTable::str(int op_code)
{
	return operator_table[op_code].str;
}
//------------------------------------------------------------------------------
inline int		OperatorTable::precedence(int op_code)
{
	return operator_table[op_code].precedence;
}
//------------------------------------------------------------------------------
inline bool	OperatorTable::isfunction(int op_code)
{
	return operator_table[op_code].isfunction;
}
//------------------------------------------------------------------------------
inline int 	OperatorTable::cost(int op_code, int what)
{
	return operator_table[op_code].cost[what];
}
//------------------------------------------------------------------------------
#endif
