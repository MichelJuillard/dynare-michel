#ifndef DYNARESCANNER_H
#define DYNARESCANNER_H
//------------------------------------------------------------------------------
/*! \file 
 \version 1.0
 \date 04/27/2004
 \par This file defines the scanner class.
*/
//------------------------------------------------------------------------------#include <iostream>
#include <sstream>
#include <string>
#include "DynareBison.h"
#ifndef YLMM_basic_scanner
#include "basic_scanner.h"
#endif
#include "Objects.h"
#include "SymbolTable.h"
#ifndef __CSTDLIB__
#include <cstdlib>
#endif
//------------------------------------------------------------------------------
/*! \namespace scanner
 */
namespace dynare 
{
  /*! 
   \class  scanner
    \brief  Member functions of this class are called from DyanreFlex.ll  
  */
  class scanner : public ylmm::basic_scanner<Objects*> 
  {
  private:
    /*! Refrence to output string */
	ostringstream				*output;
    /*! Left hand side of an expression is stored here */
	string						lhs_exp;
  public :
    /*!
     Constructor 
	 \param buf The buffer to read from  
	 */
    scanner(ylmm::basic_buffer* buf=0) 
      : ylmm::basic_scanner<Objects*>(buf) 
    {
      _current->auto_increment(true);
    }
    /*! Destructor */
    virtual ~scanner() {}
    /*! Sets reference to output string*/
    void setoutput(ostringstream* ostr)
	{
		output = ostr;	
	}
    /*! Creates a new token object name from literals string */
    void do_name(const char* name)
    {
      _token = new  Objects(name); 
    }
    /*! Creates a new token object from numerical constant*/
    void do_num_constant(const char* name)
    {
      _token = new  Objects(name); 
    }    
    /*! Creates a new token object from operator */
    void do_operator(int op)
    {
      _token = new  Objects(op); 
    }
    /*! Writes comments to output string */
    void do_comments(const char* cmt)
    {
    	*output << "% " << cmt << endl;
    }
    /*! Write token as it apears in model file */
    void do_as_is(const char* str)
    {
    	*output << str;
    }
    /*! 
     Writes symbol after removing blanks befor and after it. 
     If symbol is in symbol table it's writen as a matlab predefined variable
    */
    void do_symbol(const char* symbol)
    {
    	lhs_exp = symbol;
    	// Remove blank, newline and tab in left of symbol
    	while (lhs_exp[0] == ' ' || lhs_exp[0] =='\t'|| lhs_exp[0] =='\r' || lhs_exp[0] =='\n')
    	{
    		lhs_exp.erase(lhs_exp.begin());
    	}
    	// Remove blank, newline and tab in right of symbol
    	int s=lhs_exp.size()-1;
    	while (lhs_exp[s] == ' ' || lhs_exp[s] =='\t' || lhs_exp[s] =='\r'|| lhs_exp[s] =='\n')
    	{
    		lhs_exp.erase(lhs_exp.end()-1);
    		s=lhs_exp.size()-1;
    	}
    	int id = SymbolTable::getID(lhs_exp);
    	if (id != -1)
    	{
    		Type type = SymbolTable::getType(lhs_exp);
    		if (type == eParameter)
			{
				*output << "M_.params" << "(" << id+1 << ")";;
			}
			else if (type == eExogenous)
			{	
				*output << "oo_.exo_steady_state"<< "(" << id+1 << ")";
			}
			else if (type == eExogenousDet)
			{
				*output <<  "exedet_" << "(" << id+1 << ")";
			}
			else if (type == eEndogenous)
			{
				*output <<  "oo_.steady_state" << "(" << id+1 << ")";
			}
    	}
    	else
    	{
    		*output << lhs_exp;
    	}
	}
	/*!
	 Writes left hand side of an expression (a user variable) as equal to predefined matlab variable;
	 like :
	 	alpha = M_.params(1);
	 	x=oo_.exo_steady_state(2);
	 	...
	 */
	void do_lhs()
	{
    	int id = SymbolTable::getID(lhs_exp);
    	if (id != -1)
    	{
    		Type type = SymbolTable::getType(lhs_exp);
    		if (type == eParameter)
			{
				*output << lhs_exp << " = M_.params" << "(" << id+1 << ");\n";
			}
			else if (type == eExogenous)
			{	
				*output << lhs_exp << " = oo_.exo_steady_state"<< "(" << id+1 << ");\n";
			}
			else if (type == eExogenousDet)
			{
				*output << lhs_exp << " = exedet_" << "(" << id+1 << ");\n";
			}
			else if (type == eEndogenous)
			{
				*output << lhs_exp << " = oo_.steady_state" << "(" << id+1 << ");\n";
			}
    	}
	}
  }; 
}
//------------------------------------------------------------------------------
#endif
