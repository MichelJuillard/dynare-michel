#include "mex.h"
#include "matrix.h"

extern "C" {

  // mexFunction: Matlab Inerface point and the main application driver

void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
  {
// const char *base="base";
  mxArray* options_ = mexGetVariable("caller", "options_");

  int kOrder;
  mxArray *mxFldp = mxGetField(options_, 0, "order");
  if (mxIsNumeric(mxFldp))
    kOrder = (int) mxGetScalar(mxFldp);
  else
    kOrder = 1;
  
  mexPrintf("order: %d \n",kOrder);
  int new_order=23;
  mxArray *mxNewOrder=mxCreateDoubleScalar((double) new_order);
  mxSetField(options_,0,"order",  mxNewOrder);
  mexPutVariable("caller", "options_", options_);
  }
}




