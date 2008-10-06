/* Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */
	 
// k_order_perturbation.cpp : Defines the entry point for the DLL application.
//

#include "stdafx.h"

#include "k_order_perturbation.h"
#include "k_order_dynare.h"
#include "math.h"
#include "mex.h"

//#include "kord/approximation.h"

/* 
BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
    switch (ul_reason_for_call)
	{
		case DLL_PROCESS_ATTACH:
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
		case DLL_PROCESS_DETACH:
			break;
    }
    return TRUE;
}


// This is an example of an exported variable
K_ORDER_PERTURBATION_API int nK_order_perturbation=0;

// This is an example of an exported function.
K_ORDER_PERTURBATION_API int fnK_order_perturbation(void)
{
	return 42;
}

// This is the constructor of a class that has been exported.
// see k_order_perturbation.h for the class definition
CK_order_perturbation::CK_order_perturbation()
{ 
	return; 
}

*/

extern "C" {
	void mexFunction(int nlhs, mxArray* plhs[],
					 int nrhs, const mxArray* prhs[])
	{
		if (nrhs < 3)
			mexErrMsgTxt("Must have at least 3 input parameters.\n");
		if (nlhs == 0)
			mexErrMsgTxt("Must have at least 1 output parameter.\n");

/*		int order = (int)mxGetScalar(prhs[0]);
		if (nrhs != 12 + order) {
			mexErrMsgTxt("Must have exactly 11+order input parameters.\n");
//			return;
		}
*/
//		const mxArray* const aa = prhs[0];
/*		const mxArray* rhs1in =  (prhs[0]);
		mxArray* rhs1[1];// = rhs1in;
		rhs1[0] = 	const_cast	<mxArray* >(rhs1in);
*/
//		const mxArray* const mFname = prhs[1];
//		const char*  mFname = mxArrayToString(prhs[nrhs-1]);

	char*  mFname = mxArrayToString(prhs[nrhs-1]);


  double *y, *x, *params;
  double *residual, *g1, *g2;
  int nb_row_x, it_;

  /* Create a pointer to the input matrix y. */
  y = mxGetPr(prhs[0]);

  /* Create a pointer to the input matrix x. */
  x = mxGetPr(prhs[1]);

  /* Create a pointer to the input matrix params. */
  params = mxGetPr(prhs[2]);

  /* Fetch time index */
  it_ = (int) mxGetScalar(prhs[3]) - 1;

  /* Gets number of rows of matrix x. */
  nb_row_x = mxGetM(prhs[1]);

  residual = NULL;
  if (nlhs >= 1)
  {
     /* Set the output pointer to the output matrix residual. */
     plhs[0] = mxCreateDoubleMatrix(16,1, mxREAL);
     /* Create a C pointer to a copy of the output matrix residual. */
     residual = mxGetPr(plhs[0]);
  }

  g1 = NULL;
  if (nlhs >= 2)
  {
     /* Set the output pointer to the output matrix g1. */
     plhs[1] = mxCreateDoubleMatrix(16, 31, mxREAL);
     /* Create a C pointer to a copy of the output matrix g1. */
     g1 = mxGetPr(plhs[1]);
  }

  g2 = NULL;
 if (nlhs >= 3)
  {
     /* Set the output pointer to the output matrix g2. */
     plhs[2] = mxCreateDoubleMatrix(16, 961, mxREAL);
     /* Create a C pointer to a copy of the output matrix g1. */
     g2 = mxGetPr(plhs[2]);
  }

  /* Call the C subroutines. */



/*
		int nstat = (int)mxGetScalar(prhs[1]);
		int npred = (int)mxGetScalar(prhs[2]);
		int nboth = (int)mxGetScalar(prhs[3]);
		int nforw = (int)mxGetScalar(prhs[4]);
		int nexog = (int)mxGetScalar(prhs[5]);

		const mxArray* const ystart = prhs[6];
		const mxArray* const shocks = prhs[7];
		const mxArray* const vcov = prhs[8];
		int seed = (int)mxGetScalar(prhs[9]);
		const mxArray* const ysteady = prhs[10];
		const int* const ystart_dim = mxGetDimensions(ystart);
		const int* const shocks_dim = mxGetDimensions(shocks);
		const int* const vcov_dim = mxGetDimensions(vcov);
		const int* const ysteady_dim = mxGetDimensions(ysteady);

		int ny = nstat + npred + nboth + nforw;
		if (ny != ystart_dim[0])
			mexErrMsgTxt("ystart has wrong number of rows.\n");
		if (1 != ystart_dim[1])
			mexErrMsgTxt("ystart has wrong number of cols.\n");
		int nper = shocks_dim[1];
		if (nexog != shocks_dim[0])
			mexErrMsgTxt("shocks has a wrong number of rows.\n");
		if (nexog != vcov_dim[0])
			mexErrMsgTxt("vcov has a wrong number of rows.\n");
		if (nexog != vcov_dim[1])
			mexErrMsgTxt("vcov has a wrong number of cols.\n");
		if (ny != ysteady_dim[0])
			mexErrMsgTxt("ysteady has wrong number of rows.\n");
		if (1 != ysteady_dim[1])
			mexErrMsgTxt("ysteady has wrong number of cols.\n");
*/
//		mxArray* res = mxCreateDoubleMatrix(ny, nper, mxREAL);
/*
		try {

		// Call	int mexCallMATLAB(int nlhs, mxArray *plhs[], int nrhs,
		//			mxArray *prhs[], const char *name);
			int success = mexCallMATLAB( nlhs, plhs, nrhs-1,  rhs1 , mFname);
//			int success = mexCallMATLAB( nlhs, plhs, nrhs-1, (struct mxArray *[]) &aa , mFname);
			//plhs[0] = res;

		


		} catch (...) {
			mexErrMsgTxt("Err: Unknown error in Call MATLAB function .\n");
			mexPrintf("MexPrintf: Unknown error in Call MATLAB function %s.\n", mFname);
			return;
		}
*/
//		basic_string sFname(mFname);
	using namespace std;
		string sFname(mFname);
		string sExt("_.dll");
		mexPrintf("MexPrintf: Call exp  %d.\n", y[0]);
		double dd = exp(y[0]);
		mexPrintf("MexPrintf: exp (%d)= %d .\n", y[0],dd);

		try {
//			typedef void * (__stdcall *DynamicFn)();
			typedef void * (DynamicFn)
				(double *y, double *x, int nb_row_x, double *params, 
				int it_, double *residual, double *g1, double *g2);
			HINSTANCE DynamicHinstance;
			//DynamicFn * pDynamicFn;
			mexPrintf("MexPrintf: Call Load run DLL %s .\n", mFname);
//			DynamicHinstance=::LoadLibraryEx(strcat(mFname,"_.dll"),NULL,DONT_RESOLVE_DLL_REFERENCES);//sExt); //"_.dll");
			DynamicHinstance=::LoadLibrary(strcat(mFname,"_.dll"));//sExt); //"_.dll");
			if (DynamicHinstance==NULL)
				//return;
				throw 1;
//			typedef void * (__stdcall *DynamicFn)();
			mexPrintf("MexPrintf: 2nd Call exp  %d.\n", y[1]);
			double dd = exp(y[1]);
			mexPrintf("MexPrintf: exp (%d)= %d .\n", y[1],dd);

			mexPrintf("MexPrintf: Call GetProcAddress  %s .\n", mFname);
			DynamicFn * pDynamicFn = 
				(DynamicFn*) ::GetProcAddress(DynamicHinstance,"Dynamic");
			if (pDynamicFn == NULL)
				//return;
				throw 2;
			mexPrintf("MexPrintf: Call Dynamic  %s .\n", mFname);
//			void * objptr = (*pDynamicFn)(y,  x, nb_row_x,  params,  it_, residual,  g1,  g2);
			try{
				(*pDynamicFn)(y,  x, nb_row_x,  params,  it_, residual,  g1,  g2);
//				if (objptr == NULL)
				//return;
			}catch (...){
				DWORD dw = GetLastError();
				mexPrintf("MexPrintf: error in Call Dynamic DLL %s, %d\n", mFname, dw);
				throw 3;
			}

		} catch (int i) {
			mexPrintf("MexPrintf: error in Load and run DLL %s , %d.\n", mFname, i);
			mexErrMsgTxt("Err: An error in Load and run DLL  .\n");
			return;
		
		} catch (...) {
			mexPrintf("MexPrintf: Unknown error in Call MATLAB function %s.\n", mFname);
			mexErrMsgTxt("Err: Unknown error in Load and run DLL  .\n");
			return;
		}
	}
}