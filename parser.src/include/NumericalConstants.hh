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
#include <map>
//------------------------------------------------------------------------------
/*!
  \class  NumericalConstants
  \brief  Handles numerical constants
*/
class NumericalConstants
{
private :
  /*! Vector of numerical constants */
  std::vector<std::string> mNumericalConstants;
  //! Map matching constants to their id
  std::map<std::string, int, std::less<std::string> > numConstantsIndex;
public :
  /*! Construcor */
  NumericalConstants();
  /*! Destructor */
  ~NumericalConstants();
  /*! Adds a constant to mNumericalConstants */
  int AddConstant(std::string iConst);
  /*! Gets a constant form mNumericalConstants */
  std::string get(int iID) const;
};
//------------------------------------------------------------------------------
#endif
