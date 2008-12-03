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
// Receives call from Dynare resol and runs instead of dr1.m
// [dr,info,M_,options_,oo_] = k_ord_perturbation(M_,options_,oo_);
// receives M_,options_,oo_ and returns dr1, in addition to M_,options_,oo_
//


#include "k_ord_dynare.h"
#include "math.h"
#include "string.h"
//#include "mex.h"
//#include "k_order_perturbation.h"



#ifdef  _MSC_VER  //&&WINDOWS 


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

#endif // _MSC_VER && WINDOWS

extern "C" {
	
	// Receives call from Dynare resol and runs instead of dr1.m
	// but you need to run set_state_space first which s not in resol!!!
	// [dr,info,M_,options_,oo_] = 
	//		k_ord_perturbation(dr,check_flag,M_,options_,oo_);
	// receives dr, check_flag,, M_,options_,oo_  and 
	// returns dr1, in addition to M_,options_,oo_
	
	void mexFunction(int nlhs, mxArray* plhs[],
		int nrhs, const mxArray* prhs[])
	{
		if (nrhs < 5)
			mexErrMsgTxt("Must have at least 5 input parameters.\n");
		if (nlhs == 0)
			mexErrMsgTxt("Must have at least 1 output parameter.\n");
		
		const mxArray* dr = prhs[0];
		const int check_flag = (int)mxGetScalar(prhs[1]);
		const mxArray* M_ = prhs[2];
		const mxArray* options_=  (prhs[3]);
		const mxArray* oo_ =  (prhs[4]);
		
		mxArray* mFname =mxGetField(M_, 0, "fname");
		if(!mxIsChar(mFname)){
			mexErrMsgTxt("Input must be of type char.");
		}
		const char*  fName = mxArrayToString(mFname);

#ifdef DEBUG		
		mexPrintf("MexPrintf 2: check_flag = %d ,  fName = %s .\n", check_flag,fName);
#endif		
        int kOrder;
		mxArray* mxFldp = mxGetField(options_, 0,"order" );
		if (mxIsNumeric(mxFldp))
			kOrder = (int)mxGetScalar(mxFldp);
		else
			kOrder = 1;
		if (nrhs != 12 + kOrder) {
			mexErrMsgTxt("Must have exactly 11+order input parameters.\n");
			return;
		}
		mxFldp 	= mxGetField(M_, 0,"params" );
		double * dparams = (double *) mxGetData(mxFldp);
		int npar = (int)mxGetN(mxFldp);
		Vector * modParams =  new Vector(dparams, npar);

		mxFldp 	= mxGetField(M_, 0,"Sigma_e" );
		dparams = (double *) mxGetData(mxFldp);
		npar = (int)mxGetN(mxFldp);
		TwoDMatrix * vCov =  new TwoDMatrix(npar, npar, dparams);


		mxFldp 	= mxGetField(oo_, 0,"steady_state" ); // use in order of declaration
//		mxFldp 	= mxGetField(dr, 0,"ys" );  // and not in order of dr.order_var
		dparams = (double *) mxGetData(mxFldp);
		npar = (int)mxGetN(mxFldp);
		Vector * ySteady =  new Vector(dparams, npar);


		mxFldp = mxGetField(dr, 0,"nstatic" );
		const int nStat = (int)mxGetScalar(mxFldp);
		mxFldp = mxGetField(dr, 0,"npred" );
		const int nPred = (int)mxGetScalar(mxFldp);
		mxFldp = mxGetField(dr, 0,"nspred" );
		const int nsPred = (int)mxGetScalar(mxFldp);
		mxFldp = mxGetField(dr, 0,"nboth" );
		const int nBoth = (int)mxGetScalar(mxFldp);
		mxFldp = mxGetField(dr, 0,"nfwrd" );
		const int nForw = (int)mxGetScalar(mxFldp);
		mxFldp = mxGetField(dr, 0,"nsfwrd" );
		const int nsForw = (int)mxGetScalar(mxFldp);

		mxFldp = mxGetField(M_, 0,"exo_nbr" );
		const int nExog = (int)mxGetScalar(mxFldp);
		mxFldp = mxGetField(M_, 0,"endo_nbr" );
		const int nEndo = (int)mxGetScalar(mxFldp);
		mxFldp = mxGetField(M_, 0,"param_nbr" );
		const int nPar = (int)mxGetScalar(mxFldp);
        // it_ should be set to M_.maximum_lag
		mxFldp = mxGetField(M_, 0,"maximum_lag" );
		const int nMax_lag = (int)mxGetScalar(mxFldp);

        const int jcols = nExog+nEndo+nsPred+nsForw; // Num of Jacobian columns
        
		mxFldp 	= mxGetField(M_, 0,"endo_names" );
		const char ** endoNames = (const char **) mxGetData(mxFldp);
		const int nendo = (int)mxGetN(mxFldp);

		mxFldp 	= mxGetField(M_, 0,"exo_names" );
		const char ** exoNames = (const char **) mxGetData(mxFldp);
		const int nexo = (int)mxGetN(mxFldp);
/******
		mxFldp 	= mxGetField(M_, 0,"param_names" );
		const char ** paramNames = (char **) mxGetData(mxFldp);
		const int npar = (int)mxGetN(mxFldp);
************/
		if ((nEndo != nendo)||(nExog != nexo)) {  //(nPar != npar)
			mexErrMsgTxt("Incorrect number of input parameters.\n");
			//return;
		}

#ifdef DEBUG		
		mexPrintf("MexPrintf 2: check_flag = %d ,  fName = %s .\n", check_flag,fName);
		mexPrintf("MexPrintf 2: nEndo = %d ,  nExo = %d .\n", nEndo,nExog);
#endif		
		/* Fetch time index */
		//		int it_ = (int) mxGetScalar(prhs[3]) - 1;

		const int nSteps =0; // Dynare++ solving steps, for time being default to 0 = deterministic steady state
		const double sstol=1.e-13; //NL solver tolerance from 

		THREAD_GROUP::max_parallel_threads = 2;//params.num_threads;

        try {
			// make journal name and journal
			std::string jName(fName); //params.basename);
			jName += ".jnl";
			Journal journal(jName.c_str());
			
			//			DynamicFn * pDynamicFn = loadModelDynamicDLL (fname);
			DynamicModelDLL dynamicDLL(fName, jcols, nMax_lag, nExog);
			
			// make KordpDynare object
			KordpDynare dynare(endoNames,  nEndo, exoNames,  nExog, nPar, // paramNames,
   			   ySteady, vCov, modParams, nStat, nPred, nForw, nBoth,
			   nSteps, kOrder, journal, dynamicDLL, sstol);
    /************			
			// make list of shocks for which we will do IRFs
			vector<int> irf_list_ind;
			if (params.do_irfs_all){
                for (int i = 0; i < dynare.nexog(); i++)
                    irf_list_ind.push_back(i);
            }
    		else
            	irf_list_ind = ((const DynareNameList&)dynare.getExogNames()).selectIndices(params.irf_list);
    ****************/			
			try {
				// intiate tensor library
				tls.init(dynare.order(),
					dynare.nstat()+2*dynare.npred()+3*dynare.nboth()+
					2*dynare.nforw()+dynare.nexog());
                
				// construct main K-order approximation class
//				FistOrderApproximation app(dynare, journal, nSteps);
				Approximation app(dynare, journal, nSteps);
                // run stochastic steady 
                app.walkStochSteady();		

			} catch (const KordException& e) {
				// tell about the exception and continue
				printf("Caught (not yet fatal) Kord exception: ");
				e.print();
				JournalRecord rec(journal);
				rec << "Solution routine not finished (" << e.get_message()
					<< "), see what happens" << endrec; 
			} catch (const TLException& e) {
				mexErrMsgTxt("Caugth TL exception.");
			} catch (SylvException& e) {
				mexErrMsgTxt("Caught Sylv exception.");
            }

			/***********************            
			std::string ss_matrix_name("K_ordp");//params.prefix);
			ss_matrix_name += "_steady_states";
			
			//		ConstTwoDMatrix(app.getSS()).writeMat4(matfd, ss_matrix_name.c_str());
			
             May be needed to
			// check the approximation
			if (params.check_along_path || params.check_along_shocks
			|| params.check_on_ellipse) {
			GlobalChecker gcheck(app, THREAD_GROUP::max_parallel_threads, journal);
			if (params.check_along_shocks)
			gcheck.checkAlongShocksAndSave(matfd, params.prefix,
			params.getCheckShockPoints(),
			params.getCheckShockScale(),
			params.check_evals);
			if (params.check_on_ellipse)
			gcheck.checkOnEllipseAndSave(matfd, params.prefix,
											 params.getCheckEllipsePoints(),
											 params.getCheckEllipseScale(),
											 params.check_evals);
											 if (params.check_along_path)
											 gcheck.checkAlongSimulationAndSave(matfd, params.prefix,
											 params.getCheckPathPoints(),
											 params.check_evals);
											 }
			*****************************/
            
        // get protected derivatives from Approximation 

//          FGSContainer* rule_ders = app.GetRuleDers();
//        	FGSContainer* rule_ders_ss = app.GetRuleDersSS();
//			TwoDMatrix* gy = new TwoDMatrix(app.GetGy());
//			TwoDMatrix* gu = new TwoDMatrix(app.GetGu());
            
        // get latest ysteady 
            double * dYsteady = (dynare.getSteady().base());
            ySteady = (Vector*)(&dynare.getSteady());
            
		} catch (const KordException& e) {
			printf("Caugth Kord exception: ");
			e.print();
			return;// e.code();
		} catch (const TLException& e) {
			printf("Caugth TL exception: ");
			e.print();
			return;// 255;
		} catch (SylvException& e) {
			printf("Caught Sylv exception: ");
			e.printMessage();
			return;// 255;
		} catch (const DynareException& e) {
			printf("Caught KordpDynare exception: %s\n", e.message());
			return;// 255;
		} catch (const ogu::Exception& e) {
			printf("Caught ogu::Exception: ");
			e.print();
			return;// 255;
        }
			
			double  *dgy, *dgu, *ysteady;
			int nb_row_x;
			
			ysteady = NULL;
			if (nlhs >= 1)
			{
				/* Set the output pointer to the output matrix ysteady. */
				plhs[0] = mxCreateDoubleMatrix(nEndo,1, mxREAL);
				/* Create a C pointer to a copy of the output ysteady. */
				ysteady = mxGetPr(plhs[0]);
			}
			
			dgy = NULL;
			if (nlhs >= 2)
			{
				/* Set the output pointer to the output matrix gy. */
				plhs[1] = mxCreateDoubleMatrix(nEndo, nEndo, mxREAL);
//				plhs[1] = (double*)(gy->getData())->base();
				/* Create a C pointer to a copy of the output matrix gy. */
				dgy = mxGetPr(plhs[1]);
			}
            
			dgu = NULL;
			if (nlhs >= 3)
			{
				/* Set the output pointer to the output matrix gu. */
				plhs[2] = mxCreateDoubleMatrix(nExog, nExog, mxREAL);
//				plhs[2] = (double*)((gu->getData())->base());
				/* Create a C pointer to a copy of the output matrix gu. */
				dgu = mxGetPr(plhs[2]);
			}
            
            
			
			/* Call the C subroutines. */
//		delete params;
//		delete vCov;
//		delete ySteady;
    };
		
// Load  <model>_Dyamic DLL );
}; // end of extern C

DynamicModelDLL::DynamicModelDLL(const char * modName, const int jcols, 
								 const int nMax_lag, const int nExog)
	:jcols( jcols), nMax_lag(nMax_lag), nExog(nExog)
{
    char fName[MAX_MODEL_NAME];
    strcpy(fName,modName);
	using namespace std;
	//		string sFname(mFname);
	//		string sFname(fname);
	//		string sExt("_.dll");
//	mexPrintf("MexPrintf: Call exp  %d.\n", y[0]);
		mexPrintf("MexPrintf: Call Load run DLL %s .\n", fName);
	
	try {
		//			typedef void * (__stdcall *DynamicFn)();
		
#ifdef WINDOWS
		
		HINSTANCE dynamicHinstance;
//		dynamicHinstance=::LoadLibraryEx(strcat(fNname,"_.dll"),NULL,DONT_RESOLVE_DLL_REFERENCES);//sExt); //"_.dll");
		dynamicHinstance=::LoadLibrary(strcat(fName,"dynamic.dll"));//sExt); //"_.dll");
		if (dynamicHinstance==NULL)
			throw 1; //alt: return;
		//		(DynamicFn*)	typedef void * (__stdcall *DynamicFn)();
		mexPrintf("MexPrintf: Call GetProcAddress  %s .\n", fName);
		Dynamic = (DynamicFn *) ::GetProcAddress(dynamicHinstance,"Dynamic");
        
# else // __linux__
		
		void *dynamicHinstance = dlopen(strcat(fNname,"_dynamic.so"), RTLD_NOW);
		if((dynamicHinstance == NULL) || dlerr()){
			cerr << dlerror() << endl;
			mexPrintf("MexPrintf:Error loading DLL: %s", dlerror);
			throw 1;
		}
		void *mkr = dlsym(dynamicHinstance, "Dynamic");
		if((mkr  == NULL) || dlerr()){
			cerr << dlerror() << endl;
			mexPrintf("MexPrintf:Error finding DLL function: %s", dlerror);
			throw 2;
		}
		//The pointer to maker must be of type void *, since that is the type returned 
		//DynamicFn * 
		Dynamic = static_cast<DynamicFn*()>(mkr)();
# endif
		
//		if (Dynamic == NULL)
//			throw 3; //return;
		
	} catch (int i) {
		mexPrintf("MexPrintf: error in Load and run DLL %s , %d.\n", fName, i);
		mexErrMsgTxt("Err: An error in Load and run DLL  .\n");
		return;
		
	} catch (...) {
		mexPrintf("MexPrintf: Unknown error in Call MATLAB function %s.\n", fName);
		mexErrMsgTxt("Err: Unknown error in Load and run DLL  .\n");
		return;
	}
}

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
			


// close DLL: If the referenced object was successfully closed, 
// close() returns 0, non 0 otherwise
int DynamicModelDLL::close(){
#ifdef WINDOWS
    // MS FreeLibrary returns non 0 if OK, 0 if fails.
	bool rb=FreeLibrary(dynamicHinstance);
	if (rb)
		return 0;
	else 
		return 1;
# else // linux
	//If OK, dlclose() returns 0, non 0 otherwise
	return dlclose(dynamicHinstance);
# endif
};


void DynamicModelDLL::eval(const Vector&y, const TwoDMatrix&x, const  Vector* modParams, 
		int it_, Vector&residual, TwoDMatrix*g1, TwoDMatrix*g2){

//		DynamicDLLfunc(double *y, double *x, int nb_row_x, double *params, 
//		int it_, double *residual, double *g1, double *g2)
//        const double *dy, *dx, dbParams;
        double  *dresidual, *dg1=NULL, *dg2=NULL; 
        int length=y.length();

        if (residual.length()<length){ // dummy or insufficient
            Vector*tempv= new Vector(length );
			residual=*tempv;
            delete tempv;
		}
        if (g1!=NULL){
            if (g1->nrows()!=length){ // dummy
                delete g2;
                g1=	new TwoDMatrix( length, jcols); // and get a new one
            }
            dg1= const_cast<double*>(g1->base());
        }
        if (g2!=NULL){
            if (g2->nrows()!=length){ // dummy 
                delete g2;
                g2=	new TwoDMatrix( length, jcols*jcols);// and get a new one
            }
            dg2= const_cast<double*>(g2->base());
        }

        dresidual=const_cast<double*>(residual.base());
        double *dy=const_cast<double*>(y.base());
        double *dx=const_cast<double*>(x.base());
        double *dbParams=const_cast<double*>(modParams->base());
        
        try{
            Dynamic(dy, dx, nExog, dbParams, it_, dresidual, dg1, dg2);
        }catch (...){
            mexPrintf("MexPrintf: error in run Dynamic DLL \n");
        }
        
//		g1=	&(TwoDMatrix(double* d, int rows, int cols));
//        g1=	new TwoDMatrix( length,  sizeof(dg1)/length, dg1);
//        g2=	new TwoDMatrix( length,  sizeof(dg2)/length, dg2);
//        residual=Vector(dresidual,length );
};

void DynamicModelDLL::eval(const Vector&y, const TwoDMatrix&x, const Vector * modParams, 
		Vector&residual, TwoDMatrix*g1, TwoDMatrix*g2){

		eval(y, x, modParams, nMax_lag, residual, g1, g2);	
};

void DynamicModelDLL::eval(const Vector&y, const Vector&x, const Vector * modParams, 
		Vector&residual, TwoDMatrix*g1, TwoDMatrix*g2){

		/** ignore given exogens and create new 2D x matrix since 
		* when calling <model>_dynamic(z,x,params,it_) x must be equal to
		* zeros(M_.maximum_lag+1,M_.exo_nbr)
		**/
		const TwoDMatrix&mx = *(new const TwoDMatrix(nMax_lag+1, nExog));  

		eval(y, mx, modParams, nMax_lag, residual, g1, g2);	
};
