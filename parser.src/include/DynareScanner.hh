#ifndef _DYNARESCANNER_HH
#define _DYNARESCANNER_HH
//------------------------------------------------------------------------------
/*! \file 
 \version 1.0
 \date 04/27/2004
 \par This file defines the scanner class.
*/
//------------------------------------------------------------------------------
#include <iostream>
#include <sstream>
#include <string>
#include "DynareBison.hh"
#ifndef YLMM_basic_scanner
#include "ylmm/basic_scanner.hh"
#endif
#include "Objects.hh"
#include "SymbolTable.hh"
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
    /*! Write token as it apears in model file */
    void do_as_is(const char* str)
    {
    	*output << str;
    }
  }; 
}
//------------------------------------------------------------------------------
#endif
