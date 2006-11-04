#ifndef _NUMERICALCONSTANTS_HH
#define _NUMERICALCONSTANTS_HH
//------------------------------------------------------------------------------
/** \file
 * \version 1.0
 * \date 12/01/2003
 * \par This file defines NumericalConstants class.
 */
//------------------------------------------------------------------------------
#include <string>
#include <vector>
//------------------------------------------------------------------------------
/*!
  \class  NumericalConstants
  \brief  Handles numerical constants
*/
class NumericalConstants
{
private :
  /*! Vector of numerical constants */
  static std::vector<std::string>  mNumericalConstants;
public :
  /*! Construcor */
  NumericalConstants();
  /*! Destructor */
  ~NumericalConstants();
  /*! Adds a constant to mNumericalConstants */
  static int    AddConstant(std::string iConst);
  /*! Gets a constant form mNumericalConstants */
  static std::string  get(int iID);
};
//------------------------------------------------------------------------------
#endif
