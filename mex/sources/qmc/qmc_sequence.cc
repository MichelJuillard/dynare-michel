/*
 * Copyright (C) 2009 Dynare Team
 *
 * This file is part of Dynare (can be used outside Dynare).
 *
 * Dynare is free software: you can redistribute it and/or modify
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

#include <string.h>
#include "matrix.h"
#include "mex.h"

#if defined(__linux__)
  #define SOBOL sobol_
  #define INITSOBOL initsobol_
#else
  #define SOBOL sobol
  #define INITSOBOL initsobol
#endif

#define maxbit 30

extern "C"{
  int SOBOL(double*, unsigned int*, unsigned int*, double*, unsigned int*, unsigned int*, unsigned int*, unsigned int*);
  int INITSOBOL(unsigned int*, double*, unsigned int*, unsigned int*, unsigned int*, unsigned int*, unsigned int*);
}

void sobol_init(unsigned int* dimension, unsigned int* flag, unsigned int* seed, unsigned int* iteration, double* quasi, unsigned int* LL, unsigned int* SV)
{
  INITSOBOL(dimension, quasi, LL, iteration, SV, flag, seed);
}

void sobol_array(unsigned int* number_of_simulations, unsigned int* dimension, unsigned int* transform, 
		 unsigned int* LL, unsigned int* SV, double *sArray, unsigned int* iteration, double *quasi)
{
  SOBOL(sArray, number_of_simulations, dimension, quasi, LL, iteration, SV, transform);  
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*
  ** CASE 1 (nrhs==3 and nlhs==1) Initialization:
  **
  ** prhs[0] ==> [integer scalar]   dimension
  ** prhs[1] ==> [integer scalar]   flag {0,1,2,3}
  ** prhs[2] ==> [integer scalar]   seed 
  ** prhs[3] ==> [integer scalar]   transform {0,1}
  **
  ** plhs[0] ==> [matlab structure] qmc
  **            
  **     qmc.dimension
  **     qmc.flag
  **     qmc.transform
  **     qmc.iteration
  **     qmc.table_of_direction_numbers
  **     qmc.table_common_denominator
  **     qmc.last
  **     qmc.seed
  **
  ** CASE 2 (nrhs==2 and nlhs=2) 
  **
  ** prhs[0] ==> [matlab structure]  qmc 
  ** prhs[1] ==> [integer scalar]    sequence_size
  **
  ** plhs[0] ==> [matrix of doubles] sArray (sequence_size*dimension)     
  ** plhs[1] ==> [matlab structure]  qmc
  */

  // Get and Check input and output:
  if ( !(nrhs == 4 | nrhs == 2) )
    {
      mexErrMsgTxt("Four or two input arguments are required!");
    }
  if ( !(nlhs == 1 | nlhs == 2) )
    {
      mexErrMsgTxt("The number of output arguments should be two or one!");
    }
  if (nrhs==4)// CASE 1.
    {
      if (nlhs>1)
	{
	  mexErrMsgTxt("The number of output arguments should be one!");
	}	
      // Test of the first input argument (type):
      if (  !( mxIsNumeric(prhs[0]) & mxIsClass(prhs[0],"uint32") ) )
	{
	  mexPrintf("\t First input (dimension) has to be an integer [int32]. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      // Test of the second input argument (type):
      if (  !( mxIsNumeric(prhs[1]) & mxIsClass(prhs[1],"uint32") ) )
	{
	  mexPrintf("\t First input (flag) has to be an integer [int32]. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      // Test of the third input argument (type):
      if (  !( mxIsNumeric(prhs[1]) & mxIsClass(prhs[1],"uint32") ) )
	{
	  mexPrintf("\t First input (seed) has to be an integer [int32]. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      // Test of the fourth input argument (type):
      if (  !( mxIsNumeric(prhs[1]) & mxIsClass(prhs[1],"uint32") ) )
	{
	  mexPrintf("\t First input (transform) has to be an integer [int32]. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}

      /*
      ** Get the input variables.
      */

      unsigned int dimension = (unsigned int) mxGetScalar(prhs[0]);
      if (dimension>1111 | dimension < 2)
	{
	  mexErrMsgTxt("dimension has to be between 2 and 1111!");
	}
      mxArray* DIMENSION = mxDuplicateArray(prhs[0]);
      
      unsigned int flag = (unsigned int) mxGetScalar(prhs[1]);
      if ( !( flag==0 | flag==1 | flag==2 | flag==3 ) )
	{
	  mexErrMsgTxt("flag has to be equal to 0, 1, 2 or 3!");
	}
      mxArray* FLAG = mxDuplicateArray(prhs[1]);
      
      unsigned int seed = ( unsigned int) mxGetScalar( prhs[2] );

      unsigned int transform = (unsigned int) mxGetScalar(prhs[3]);
      if ( !( transform==0 | transform==1) )
	{
	  mexErrMsgTxt("transform has to be equal to 0 or 1!");
	}
      mxArray* TRANSFORM = mxDuplicateArray(prhs[3]);

      /*
      ** Declaration of the variables returned by INITSOBOL.
      */

      mxArray* t0;
      unsigned int* LL = (unsigned int*) mxMalloc(sizeof(mxUINT32_CLASS));
      t0 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
      LL = (unsigned int*) mxGetData(t0);

      mxArray* t1;
      unsigned int* SV = (unsigned int*) mxMalloc(dimension*maxbit*sizeof(mxUINT32_CLASS));
      t1 = mxCreateNumericMatrix(dimension, maxbit, mxUINT32_CLASS, mxREAL);
      SV = (unsigned int*) mxGetData(t1);

      mxArray* t2;
      unsigned int* iteration = (unsigned int*) mxMalloc((sizeof(mxUINT32_CLASS))); 
      t2 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
      iteration = (unsigned int*) mxGetData(t2);
      mxArray* t3;

      double* quasi = (double*) mxMalloc((dimension*sizeof(double)));
      t3 = mxCreateNumericMatrix(dimension,1,mxDOUBLE_CLASS,mxREAL);
      quasi = (double*) mxGetData(t3);

      mxArray* t4;
      unsigned int* seedout = (unsigned int*) mxMalloc(sizeof(mxUINT32_CLASS));
      t4 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS,mxREAL);
      seedout = (unsigned int*) mxGetData(t4);

      /*
      ** Call to INITSOBOL subroutine.
      */

      sobol_init(&dimension, &flag, &seed, iteration, quasi, LL, SV);

      mxSetData(t0,LL);
      mxSetData(t1,SV);
      mxSetData(t2,iteration);
      mxSetPr(t3,quasi);
      *seedout = seed;
      mxSetData(t4,seedout);

      /*
      ** Now I build the returned matlab structure.
      */

      char *fieldnames[8]; //This will hold field names.
      fieldnames[0] = (char*) mxMalloc(sizeof("dimension"));
      memcpy(fieldnames[0],"dimension",sizeof("dimension"));
      fieldnames[1] = (char*) mxMalloc(sizeof("flag"));
      memcpy(fieldnames[1],"flag",sizeof("flag"));
      fieldnames[2] = (char*) mxMalloc(sizeof("transform"));
      memcpy(fieldnames[2],"transform",sizeof("transform"));
      fieldnames[3] = (char*) mxMalloc(sizeof("iteration"));
      memcpy(fieldnames[3],"iteration",sizeof("iteration"));
      fieldnames[4] = (char*) mxMalloc(sizeof("table_of_direction_numbers"));
      memcpy(fieldnames[4],"table_of_direction_numbers",sizeof("table_of_direction_numbers"));
      fieldnames[5] = (char*) mxMalloc(sizeof("table_common_denominator"));
      memcpy(fieldnames[5],"table_common_denominator",sizeof("table_common_denominator"));
      fieldnames[6] = (char*) mxMalloc(sizeof("last"));
      memcpy(fieldnames[6],"last",sizeof("last"));
      fieldnames[7] = (char*) mxMalloc(sizeof("seed"));
      memcpy(fieldnames[7],"seed",sizeof("seed"));

      plhs[0] = mxCreateStructMatrix(1,1,8,(const char**)fieldnames);

      mxFree( fieldnames[0] );
      mxFree( fieldnames[1] );
      mxFree( fieldnames[2] );
      mxFree( fieldnames[3] );
      mxFree( fieldnames[4] );
      mxFree( fieldnames[5] );
      mxFree( fieldnames[6] );
      mxFree( fieldnames[7] );

      mxSetField(plhs[0], 0, "dimension", DIMENSION);
      mxSetField(plhs[0], 0, "flag", FLAG);
      mxSetField(plhs[0], 0, "transform", TRANSFORM);
      mxSetField(plhs[0], 0, "iteration", t2);
      mxSetField(plhs[0], 0, "table_of_direction_numbers", t1);
      mxSetField(plhs[0], 0, "table_common_denominator", t0);
      mxSetField(plhs[0], 0, "last", t3);
      mxSetField(plhs[0], 0, "seed", t4);
  
    }
  if (nrhs==2)// CASE 2.
    {
      if (nlhs!=2)
	{
	  mexErrMsgTxt("The number of output arguments has to be two!");
	}
      if ( !mxIsStruct(prhs[0]) )
	{
	  mexErrMsgTxt("The first argument must be a matlab structure!");
	}
      else
	{
	  if (mxGetNumberOfFields(prhs[0]) != 8)
	    {
	      mexErrMsgTxt("The matlab structure must have eight fields!");
	    }
	}

      /*
      ** Get the fields of the matlab structure.
      */

      plhs[1] = mxDuplicateArray(prhs[0]);

      mxArray* t0;// dimension
      t0 = mxGetField( plhs[1], 0, "dimension");
      if (t0 == NULL)
	{
	  mexPrintf("\t Field ""dimension"" is empty. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      if (!mxIsClass(t0,"uint32"))
	{
	  mexPrintf("\t Field ""dimension"" is not uint32. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      unsigned int* dimension = (unsigned int*) mxMalloc(sizeof( mxUINT32_CLASS ));
      dimension = (unsigned int*) mxGetData( t0 );

//       mxArray* t1;// flag
//       t1 = mxGetField(plhs[1], 0, "flag");
//       if (t1 == NULL)
// 	{
// 	  mexPrintf("\t Field ""flag"" is empty. \n");
// 	  mexErrMsgTxt("\t Fatal error.");
// 	}
//       if (!mxIsClass(t1,"uint32"))
// 	{
// 	  mexPrintf("\t Field ""flag"" is not uint32. \n");
// 	  mexErrMsgTxt("\t Fatal error.");
// 	}
//       unsigned int* flag = (unsigned int*) mxMalloc(sizeof( mxUINT32_CLASS ));
//       flag = (unsigned int*) mxGetData( t1 );
      
      mxArray* t2;// transform
      t2 = mxGetField(plhs[1], 0, "transform");
      if (t2 == NULL)
	{
	  mexPrintf("\t Field ""transform"" is empty. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      if (!mxIsClass(t2,"uint32"))
	{
	  mexPrintf("\t Field ""transform"" is not uint32. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      unsigned int* transform = (unsigned int*) mxMalloc(sizeof( mxUINT32_CLASS ));
      transform = (unsigned int*) mxGetData( t2 );
      
      mxArray* t3;// iteration
      t3 = mxGetField(plhs[1], 0, "iteration");
      if (t3 == NULL)
	{
	  mexPrintf("\t Field ""iteration"" is empty. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      if (!mxIsClass(t3,"uint32"))
	{
	  mexPrintf("\t Field ""iteration"" is not uint32. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      unsigned int* iteration = (unsigned int*) mxMalloc(sizeof( mxUINT32_CLASS ));
      iteration = (unsigned int*) mxGetData( t3 );
      
      mxArray* t4;// table_of_direction_numbers
      t4 = mxGetField(plhs[1], 0, "table_of_direction_numbers");
      if (t4 == NULL)
	{
	  mexPrintf("\t Field ""table_of_direction_numbers"" is empty. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      if (!mxIsClass(t4,"uint32"))
	{
	  mexPrintf("\t Field ""table_of_direction_numbers"" is not uint32. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      unsigned int* SV = (unsigned int*) mxMalloc(*dimension*maxbit*sizeof( mxUINT32_CLASS ));
      SV = (unsigned int*) mxGetData( t4 );
      
      mxArray* t5;//table_common_denominator
      t5 = mxGetField(plhs[1], 0, "table_common_denominator");
      if (t5 == NULL)
	{
	  mexPrintf("\t Field ""table_common_denominator"" is empty. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      if (!mxIsClass(t5,"uint32"))
	{
	  mexPrintf("\t Field ""table_common_denominator"" is not uint32. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      unsigned int* LL = (unsigned int*) mxMalloc(sizeof( mxUINT32_CLASS ));
      LL = (unsigned int*) mxGetData( t5 );
      
      mxArray* t6;// last
      t6 = mxGetField(plhs[1], 0, "last");
      if (t6 == NULL)
	{
	  mexPrintf("\t Field ""last"" is empty. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      if (!mxIsClass(t6,"double"))
	{
	  mexPrintf("\t Field ""last"" is not double. \n");
	  mexErrMsgTxt("\t Fatal error.");
	}
      double* last = (double*) mxMalloc(*dimension*sizeof(double));
      last = (double*) mxGetData( t6 );
      
//       mxArray* t7;// seed
//       t7 = mxGetField(prhs[0], 0, "seed");
//       if (t7 == NULL)
// 	{
// 	  mexPrintf("\t Field ""seed"" is empty. \n");
// 	  mexErrMsgTxt("\t Fatal error.");
// 	}
//       if (!mxIsClass(t7,"uint32"))
// 	{
// 	  mexPrintf("\t Field ""seed"" is not uint32. \n");
// 	  mexErrMsgTxt("\t Fatal error.");
// 	}
//       unsigned int* seed = (unsigned int*) mxMalloc(sizeof( mxUINT32_CLASS ));
//       seed = (unsigned int*) mxGetData( t7 );
      
      /*
      ** Get the second input.
      */
      
      unsigned int number_of_simulations = (unsigned int) mxGetScalar(prhs[1]); 


      /*
      ** Initialization of the first output argument.
      */

      plhs[0] = mxCreateNumericMatrix(number_of_simulations, *dimension, mxDOUBLE_CLASS, mxREAL);
      double* sArray = (double*)mxMalloc(*dimension*number_of_simulations*sizeof(double));
      mxSetPr(plhs[0],sArray);

      /*
      ** Call to SOBOL subroutine.
      */

      sobol_array(&number_of_simulations, dimension, transform, LL, SV, sArray, iteration, last);

      mxSetData(t6,last);
      mxSetField(plhs[1], 0, "last", t6);

      mxSetData(t3,iteration);
      mxSetField(plhs[1], 0, "iteration", t3);

     }
}
