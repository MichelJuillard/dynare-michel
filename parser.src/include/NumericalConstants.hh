#ifndef _NUMERICALCONSTANTS_HH
#define _NUMERICALCONSTANTS_HH

using namespace std;

#include <string>
#include <vector>
#include <map>

//! Handles numerical constants
class NumericalConstants
{
private:
  //! Vector of numerical constants
  vector<string> mNumericalConstants;
  //! Map matching constants to their id
  map<string, int> numConstantsIndex;
public:
  NumericalConstants();
  //! Adds a constant and returns its ID
  int AddConstant(const string &iConst);
  //! Get a constant in string form
  string get(int iID) const;
  //! Get a constant in double form
  double getDouble(int iID) const;
};

#endif
