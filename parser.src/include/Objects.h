#ifndef OBJECTS_HH
#define OBJECTS_HH
//------------------------------------------------------------------------------
/*! \file 
 \version 1.0
 \date 04/13/2004
 \par This file defines the Objects class.
*/
//------------------------------------------------------------------------------
using namespace std;
#include "SymbolTable.h"
#include "ModelTypes.h"
//------------------------------------------------------------------------------
/*!
  \class  Objects
  \brief  This class defines data associeted to parsed tokens.
*/
namespace dynare 
{
	class Objects
	{
		public :
			/*! Parsed string name : can be literal name or operator */
			string 			symbol;
			/*! ID of object : can be ID of symbol, variable, or token in model tree */
			NodeID 			ID;
			/*! Type of object */
			Type			type;			
			/*! In case of operator object, this is set to its code */
			int				opcode;
		public :
			/*! Constrcutor with default values */
			Objects() 
			{
				ID = NULL;				
				symbol = "";
				type = eUNDEF;
			}
			/*! Constrcutor of object with known symbol name */
			Objects(string name) 
			{
				symbol = name;					
			}
			/*! Constrcutor of object with known symbol name */
			Objects(const char* name, NodeID id = NULL, Type t = eUNDEF)  
			{
				symbol = name;					
				type = t;	
				ID = id;
			}
			/*! Constrcutor of object with known symbol name, ID and type */
			Objects(string name,NodeID id, Type t) 
			{
				symbol = name;
				type = t;	
				ID = id;					
			}
			/*! Constrcutor of object with known operator code */
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
