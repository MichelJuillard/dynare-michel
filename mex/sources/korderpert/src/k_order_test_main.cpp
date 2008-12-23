// k_order_test_main.cpp 

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

/*************************************
* This main() is for testing k_order DLL entry point by linking to 
* the k_ord library statically and passing its hard-coded data: 
* parameters, covar, ysteady and the variable names from fs2000a.mod model
* The main has been derived from mxFunction used for K-Order DLL
***************************************/

//#include "stdafx.h"
#include "k_ord_dynare.h"


int main(int argc, char* argv[])
{


	const int check_flag = 0;
	const char*  fName = "fs2000a";//mxArrayToString(mFname);
	
#ifdef DEBUG		
	mexPrintf("k_order_perturbation: check_flag = %d ,  fName = %s .\n", check_flag,fName);
#endif		
	int kOrder =1;
	int npar = 7;//(int)mxGetM(mxFldp);
	double dparams[7]={ 0.3300,
		0.9900,
		0.0030,
		1.0110,
		0.7000,
		0.7870,
		0.0200
	};
	
	Vector * modParams =  new Vector(dparams, npar);
	
#ifdef DEBUG		
    mexPrintf("k_ord_perturbation: nParams=%d .\n",npar);  
    for (int i = 0; i < npar; i++) {
        mexPrintf("k_ord_perturbation: dParams[%d]= %g.\n", i, dparams+i*(sizeof(double)) );  }
    for (int i = 0; i < npar; i++) {
        mexPrintf("k_ord_perturbation: Params[%d]= %g.\n", i, (*modParams)[i]);  }
	
#endif		        
	
	double d2Dparams[4] = { //(double *) mxGetData(mxFldp);
			0.1960e-3, 0.0,
			0.0, 0.0250e-3};
		npar = 2;//(int)mxGetN(mxFldp);
		TwoDMatrix * vCov =  new TwoDMatrix(npar, npar, (d2Dparams));
		
		double dYSparams [29]= { // mxGetData(mxFldp);
			1.0110, 2.2582, 5.8012, 1.0000, 1.0000, 0.5808, 1.0110, 2.2582,
			0.4477, 1.0000, 4.5959, 1.0212, 5.8012, 0.8494, 0.1872, 0.8604,
			1.0030, 1.0080, 1.0000, 1.0000, 0.5808, 1.0030, 1.0110, 2.2582,
		    0.4477, 1.0000, 0.1872, 2.2582, 0.4477
		};
//			1.0110, 2.2582, 0.4477, 1.0000, 4.5959, 1.0212, 5.8012, 0.8494,
//				0.1872, 0.8604, 1.0030, 1.0080, 1.0000, 1.0000, 0.5808, 1.0030
//		};
		const int nSteady = 29;//16 (int)mxGetM(mxFldp);
		Vector * ySteady =  new Vector(dYSparams, nSteady);
		
		//mxFldp = mxGetField(dr, 0,"nstatic" );
		const int nStat = 7;//(int)mxGetScalar(mxFldp);
		//	mxFldp = mxGetField(dr, 0,"npred" );
		const int nPred = 6;//(int)mxGetScalar(mxFldp);
		//mxFldp = mxGetField(dr, 0,"nspred" );
		const int nsPred = 6;//(int)mxGetScalar(mxFldp);
		//mxFldp = mxGetField(dr, 0,"nboth" );
		const int nBoth = 2;// (int)mxGetScalar(mxFldp);
		//mxFldp = mxGetField(dr, 0,"nfwrd" );
		const int nForw = 3;// (int)mxGetScalar(mxFldp);
		//mxFldp = mxGetField(dr, 0,"nsfwrd" );
		const int nsForw = 7;//(int)mxGetScalar(mxFldp);
		
		//mxFldp = mxGetField(M_, 0,"exo_nbr" );
		const int nExog =2;// (int)mxGetScalar(mxFldp);
		//mxFldp = mxGetField(M_, 0,"endo_nbr" );
		const int nEndo = 16;//(int)mxGetScalar(mxFldp);
		//mxFldp = mxGetField(M_, 0,"param_nbr" );
		const int nPar = 7;//(int)mxGetScalar(mxFldp);
        // it_ should be set to M_.maximum_lag
		//mxFldp = mxGetField(M_, 0,"maximum_lag" );
		const int nMax_lag = 1;//(int)mxGetScalar(mxFldp);
		
        const int jcols = nExog+nEndo+nsPred+nsForw; // Num of Jacobian columns
        mexPrintf("k_order_perturbation: jcols= %d .\n", jcols);
		
        //mxFldp= mxGetField(M_, 0,"endo_names" );
        const int nendo = 16;//(int)mxGetM(mxFldp);
        const int widthEndo = 6;// (int)mxGetN(mxFldp);
		const char * cNamesCharStr= "mPceWRkdnlggYPyd          yp__ A          __oo            oobb            bbss            ss    ";
		//       const char**  endoNamesMX= DynareMxArrayToString( mxFldp,nendo,widthEndo);
		const char**  endoNamesMX= DynareMxArrayToString( cNamesCharStr, nendo, widthEndo);
#ifdef DEBUG		
        for (int i = 0; i < nEndo; i++) {
            mexPrintf("k_ord_perturbation: EndoNameList[%d][0]= %s.\n", i, endoNamesMX[i] );
        }
#endif	     
        //mxFldp 	= mxGetField(M_, 0,"exo_names" );
		const int nexo = 2;//(int)mxGetM(mxFldp);
		const int widthExog = 3;//(int)mxGetN(mxFldp);
		//        const char**  exoNamesMX= DynareMxArrayToString( mxFldp,nexo,widthExog);
		const char * cExoNamesCharStr= "ee__am";
        const char**  exoNamesMX= DynareMxArrayToString( cExoNamesCharStr,nexo,widthExog);
        
#ifdef DEBUG		
        for (int i = 0; i < nexo; i++) {
            mexPrintf("k_ord_perturbation: ExoNameList[%d][0]= %s.\n", i, exoNamesMX[i] );
        }
#endif	     
		if ((nEndo != nendo)||(nExog != nexo)) {  //(nPar != npar)
			mexErrMsgTxt("Incorrect number of input parameters.\n");
			//return;
		}
		
#ifdef DEBUG		
		for (int i = 0; i < nEndo; i++) {
			mexPrintf("k_ord_perturbation: EndoNameList[%d]= %s.\n", i, endoNamesMX[i] );   }
		//	for (int i = 0; i < nPar; i++) {
		//        mexPrintf("k_ord_perturbation: params_vec[%d]= %g.\n", i, params_vec[i] );   }
		for (int i = 0; i < nPar; i++) {
			mexPrintf("k_ord_perturbation: Params[%d]= %g.\n", i, (*modParams)[i]);  }
		for (int i = 0; i < nSteady; i++) {
			mexPrintf("k_ord_perturbation: ysteady[%d]= %g.\n", i, (*ySteady)[i]);  }
		
		mexPrintf("k_order_perturbation: nEndo = %d ,  nExo = %d .\n", nEndo,nExog);
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
#ifdef DEBUG		
			mexPrintf("k_order_perturbation: Calling dynamicDLL constructor.\n");
#endif				
			//			DynamicFn * pDynamicFn = loadModelDynamicDLL (fname);
			DynamicModelDLL dynamicDLL(fName, nEndo, jcols, nMax_lag, nExog);
#ifdef DEBUG		
			mexPrintf("k_order_perturbation: Calling dynare constructor.\n");
#endif			
			// make KordpDynare object
			KordpDynare dynare(endoNamesMX,  nEndo, exoNamesMX,  nExog, nPar, // paramNames,
				ySteady, vCov, modParams, nStat, nPred, nForw, nBoth,
				nSteps, kOrder, journal, dynamicDLL, sstol);
			try {
				// intiate tensor library
#ifdef DEBUG		
				mexPrintf("k_order_perturbation: Call tls init\n");
#endif
                tls.init(dynare.order(),
					dynare.nstat()+2*dynare.npred()+3*dynare.nboth()+
					2*dynare.nforw()+dynare.nexog());
                
				// construct main K-order approximation class
				//				FistOrderApproximation app(dynare, journal, nSteps);
#ifdef DEBUG		
				mexPrintf("k_order_perturbation: Call Approximation constructor \n");
#endif
				Approximation app(dynare, journal, nSteps);
                // run stochastic steady 
#ifdef DEBUG		
				mexPrintf("k_order_perturbation: Calling walkStochSteady.\n");
#endif			
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
			
            
			// get latest ysteady 
            double * dYsteady = (dynare.getSteady().base());
            ySteady = (Vector*)(&dynare.getSteady());
            
		} catch (const KordException& e) {
			printf("Caugth Kord exception: ");
			e.print();
			return 1;// e.code();
		} catch (const TLException& e) {
			printf("Caugth TL exception: ");
			e.print();
			return 2;// 255;
		} catch (SylvException& e) {
			printf("Caught Sylv exception: ");
			e.printMessage();
			return 3;// 255;
		} catch (const DynareException& e) {
			printf("Caught KordpDynare exception: %s\n", e.message());
			return 4;// 255;
		} catch (const ogu::Exception& e) {
			printf("Caught ogu::Exception: ");
			e.print();
			return 5;// 255;
		}

		// bones for future developement of the output.

		const int nrhs=5;
		const int nlhs=2;

		mxArray* prhs[nrhs];
		mxArray* plhs[nlhs];

#ifdef DEBUG		
		mexPrintf("k_order_perturbation: Filling outputs.\n");
#endif			
		
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
			plhs[1] = mxCreateDoubleMatrix(nEndo, jcols, mxREAL);
			//				plhs[1] = (double*)(gy->getData())->base();
			/* Create a C pointer to a copy of the output matrix gy. */
			dgy = mxGetPr(plhs[1]);
		}
		
		dgu = NULL;
		if (nlhs >= 3)
		{
			/* Set the output pointer to the output matrix gu. */
			plhs[2] = mxCreateDoubleMatrix(nEndo, nExog, mxREAL);
			//				plhs[2] = (double*)((gu->getData())->base());
			/* Create a C pointer to a copy of the output matrix gu. */
			dgu = mxGetPr(plhs[2]);
		}
		
		
		return 0;
}
