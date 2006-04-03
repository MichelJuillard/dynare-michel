#ifndef NUMERICALINITIALIZATION_H
#define NUMERICALINITIALIZATION_H
//------------------------------------------------------------------------------
/*! \file 
 \version 1.0
 \date 04/16/2004
 \par This file defines the NumericalInitialization class.
*/
//------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <sstream>
//------------------------------------------------------------------------------
#include "SymbolTable.h"
//------------------------------------------------------------------------------
/*! \class  NumericalInitialization
 \brief  Handles numerical initialization of Endogenousous and Exogenousous variables.
*/
class NumericalInitialization
{
	private :
		/*! Output of this class */
		std::ostringstream	*output;			
	public :
		/*! Constructor */
		NumericalInitialization();
		/*! Destrcutor */
		~NumericalInitialization();
		/*! Pointer to error function of parser class */                               
		void (* error) (const char* m);
		/*! Sets reference to external output string */
		void setOutput(std::ostringstream* iOutput); 
		/*! 
		  \param name a string.
		  \param expression an Expression type.
		  \par Description
		  Set constant value
		  - in interpretative languages (Matlab, Gauss, Scilab)\n
		  - in C++, evaluate expression and set value for Name in Parameters_Table\n
		 */
		void SetConstant(std::string name,std::string expression);
		void SetLocalConstant(std::string name,std::string expression);
		/*!
		 \par Description
		 Initializes an initval block to set initial values for variables 
		 */
		void BeginInitval(void);
		/*!
		  \param name a string.
		  \param expression an Expression type.
		  \par Description
		  Set initial or terminal value for a variable : 
		  - name must be searched for in SymbolTable 
		  - error if symbol isn't a variable (Endogenousous or Exogenousous or Exogenousous deterministic or recursive) 
		  - write "oo_.steady_state(symbol.id) = expression;" (for Endogenousous variables) or "oo_.exo_steady_state(symbol.id) = expression;" (for Exogenousous variables) 
		 */ 
		void SetInit(std::string name,std::string expression);
		/*!
		 \par Description
		 Closes an initval block.
		 */
		void EndInitval(void);
		/*!
		 \par Description
		 Initializes an endval block to set terminal values for variables.
		 */
		void BeginEndval(void);
		/*!
		 \par Description
		 Closes an endval block.
		 */		
		void EndEndval(void);
		/*!
		 \par Description
		 Initializes an histval block to set historical values for variables.
		 */
		void BeginHistval(void);		
		/*!
		  \param name a string.
		  \param lag a string.
		  \param expression an Expression type.
		  \par Description
		  Set historical for a variable :
		  - name must be searched for in SymbolTable 
		  - error if symbol isn't a variable (Endogenousous or Exogenousous or Exogenousous deterministic or recursive) 
		  - write "oo_.y_simul(symbol.id,MP.max_lag+lag+1) = expression;" (for Endogenousous variables) or "oo_.exo_simul(MP.max_lag+lag+1,symbol.id) = expression;" (for Exogenousous variables) (lag is <= 0, MP.max_lag will be set in ModelTree) 
		 */
		void SetHist(std::string name,int lag, std::string expression);
};
//------------------------------------------------------------------------------
#endif
