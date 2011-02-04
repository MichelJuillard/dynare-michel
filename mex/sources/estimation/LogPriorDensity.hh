///////////////////////////////////////////////////////////
//  LogPriorDensity.hh
//  Implementation of the Class LogPriorDensity
//  Created on:      10-Feb-2010 20:56:08
///////////////////////////////////////////////////////////

#if !defined(LPD_011FD4CF_17CE_4805_882B_046AA07CF443__INCLUDED_)
#define LPD_011FD4CF_17CE_4805_882B_046AA07CF443__INCLUDED_

//#include <boost/random/variate_generator.hpp>
#include "Vector.hh"
#include "EstimatedParametersDescription.hh"

class LogPriorDensity
{

public:
  LogPriorDensity(EstimatedParametersDescription &estParsDesc);
  virtual ~LogPriorDensity();

  double compute(const Vector &estParams);
  void computeNewParams(Vector &newParams);

private:
  const EstimatedParametersDescription &estParsDesc;

};

#endif // !defined(011FD4CF_17CE_4805_882B_046AA07CF443__INCLUDED_)
