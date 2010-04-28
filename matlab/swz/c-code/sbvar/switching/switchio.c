
#include "switchio.h"
#include "dw_array.h"
#include "dw_matrix_array.h"
#include "dw_error.h"
#include "dw_ascii.h"

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "modify_for_mex.h"

static void ReadError(char *idformat, char *trailer, int error);
static int SetFilePosition(FILE *f_in, char *format, char *str);
static int ReadInteger(FILE *f_in, char *idformat, char *trailer, int *i);
static int ReadMatrix(FILE *f_in, char *idformat, char *trailer, TMatrix X);
static int ReadIntArray(FILE *f_in, char *idformat, char *trailer, void *X);

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
#define SWITCHIO_LINE_ID_NOT_FOUND  1
#define SWITCHIO_ERROR_READING_DATA 2

static void ReadError(char *idformat, char *trailer, int error)
{
  char *idbuffer, *errmsg, *format;
  switch (error)
    {
    case SWITCHIO_LINE_ID_NOT_FOUND:
      format="Line identifier ""%s"" not found.";
      break;
    case SWITCHIO_ERROR_READING_DATA:
      format="Error reading data after line identifier ""%s"".";
      break;
    case 0:
      return;
    default:
      dw_Error(UNKNOWN_ERR);
      return;
    }
  if (trailer)
    sprintf(idbuffer=(char*)malloc(strlen(idformat)+strlen(trailer)-1),idformat,trailer);
  else
    idbuffer=idformat;
  sprintf(errmsg=(char*)malloc(strlen(format)+strlen(idbuffer)-1),format,idbuffer);
  dw_UserError(errmsg);
  free(errmsg);
  if (idbuffer != idformat) free(idbuffer);
}

/*
   Assumes
     format : "*%s*" or "*" 
     str    : "*" 

     where * is any string that does not contain format specifiers.

   Results
     finds given string in file
   
*/
static int SetFilePosition(FILE *f_in, char *format, char *str)
{
  char *buffer;
  int rtrn;
  if (str)
    sprintf(buffer=(char*)malloc(strlen(format)+strlen(str)-1),format,str);
  else
    buffer=format;
  rtrn=dw_SetFilePosition(f_in,buffer);
  if (buffer != format) free(buffer);
  return rtrn;
}

static int ReadInteger(FILE *f_in, char *format, char *str, int *i)
{
  if (!SetFilePosition(f_in,format,str))
    return SWITCHIO_LINE_ID_NOT_FOUND;
  else
    if (fscanf(f_in," %d ",i) != 1)
      return SWITCHIO_ERROR_READING_DATA;
    else
      return 0;
}

static int ReadMatrix(FILE *f_in, char *format, char *str, TMatrix X)
{
  if (!SetFilePosition(f_in,format,str))
    return SWITCHIO_LINE_ID_NOT_FOUND;
  else
    if (!dw_ReadMatrix(f_in,X))
      return SWITCHIO_ERROR_READING_DATA;
    else
      return 0;
}

static int ReadIntArray(FILE *f_in, char *format, char *str, void *X)
{
  if (!SetFilePosition(f_in,format,str))
    return SWITCHIO_LINE_ID_NOT_FOUND;
  else
    if (!dw_ReadArray(f_in,X))
      return SWITCHIO_ERROR_READING_DATA;
    else
      return 0;
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/***************** TMarkovStateVariable Input/Output Routines ******************/
/*******************************************************************************/
/*
   Assumes:
    f:  valid file pointer
    idstring:  pointer to null terminated string

   Returns:
    Valid pointer to a TMarkovStateVariable structure upon success and null 
    pointer upon failure.

   Results:
    Reads Markov specification from file and creates TMarkovStateVariable 
    structure. 

   Notes:
    The null terminated string idstring is of the form "" or "[i][j]...[k]".  
    Usually, this function will be called with idstring equal to "".  The routine 
    is called recursively with the more general form for idstring.  For each path 
    of Markov state variables of the form

        sv->state_variable[i]->state_variable[j]-> ... ->state_variable[k]

    which corresponds to a single Markov state variable, there must be an entries 
    in the file of the form:
*/
TMarkovStateVariable* ReadMarkovSpecification_SV(FILE *f_in, char *idstring, int nobs)
{
  char *idformat, *trailer=(char*)NULL, *idstring_new;
  int err, i, j, n_state_variables, nstates, nlags_encoded, nbasestates;
  TMatrix Prior, MQ;
  int *FreeDim, **NonZeroIndex;
  TMarkovStateVariable *sv=(TMarkovStateVariable*)NULL, **sv_array;

  // Get nobs if necessary
  if ((nobs > 0) || !(err=ReadInteger(f_in,idformat="//== Number observations ==//",trailer,&nobs)))
    {
      // Construct trailer
      if (idstring[0])
	sprintf(trailer=(char*)malloc(24+strlen(idstring)),"for state_variable%s ==//",idstring);
      else
	strcpy(trailer=(char*)malloc(5),"==//");

      // Read number of state variables
      if (!(err=ReadInteger(f_in,idformat="//== Number independent state variables %s",trailer,&n_state_variables)))
	{
	  if (n_state_variables > 1)
	    {
	      sv_array=(TMarkovStateVariable**)dw_CreateArray_pointer(n_state_variables,(void (*)(void*))FreeMarkovStateVariable);
	      for (j=10, i=1; n_state_variables/j > 0; j*=10, i++);
	      strcpy(idstring_new=(char*)malloc((j=(int)strlen(idstring))+i+3),idstring); 
	      for (i=0; i < n_state_variables; i++)
		{
		  sprintf(idstring_new+j,"[%d]",i+1);
		  if (!(sv_array[i]=ReadMarkovSpecification_SV(f_in,idstring_new,nobs))) break;
		}
	      free(idstring_new);
	      if (i == n_state_variables)
		sv=CreateMarkovStateVariable_Multiple(nobs,n_state_variables,sv_array);
	      else
		dw_FreeArray(sv_array);
	    }
	  else
	    {
	      // Read number states
	      if (!(err=ReadInteger(f_in,idformat="//== Number states %s",trailer,&nstates)))
		{
		  // Read number of lags to encode
		  switch (err=ReadInteger(f_in,idformat="//== Number of lags encoded %s",trailer,&nlags_encoded))
		    {
		    case SWITCHIO_ERROR_READING_DATA:
		      break;
		    case SWITCHIO_LINE_ID_NOT_FOUND:
		      nlags_encoded=0;
		    case 0:
		      if (nlags_encoded < 0)
			{
			  err=SWITCHIO_ERROR_READING_DATA;
			  break;
			}

		      // Read number of base states
		      if (nlags_encoded > 0)
			{
			  if (err=ReadInteger(f_in,idformat="//== Number of base states %s",trailer,&nbasestates))
			    break;
			  for (j=nbasestates, i=nlags_encoded; i > 0; i--) j*=nbasestates;
			  if (j != nstates)
			    {
			      err=SWITCHIO_ERROR_READING_DATA;
			      break;
			    }
			}

		      // Read prior
		      Prior=CreateMatrix(nstates,nstates);
		      if (!(err=ReadMatrix(f_in,idformat="//== Prior %s",trailer,Prior)))
			{
			  // Read free Dirichlet dimensions
			  if (!(err=ReadInteger(f_in,idformat="//== Number free Dirichlet variables %s",trailer,&i)))
			    {
			      FreeDim=dw_CreateArray_int(i);
			      if (!(err=ReadIntArray(f_in,idformat="//== Free Dirichlet dimensions %s",trailer,FreeDim)))
				{
				  // Read free Dirichlet index
				  NonZeroIndex=dw_CreateRectangularArray_int(nstates,nstates);
				  if (!(err=ReadIntArray(f_in,idformat="//== Free Dirichlet index %s",trailer,NonZeroIndex)))
				    {
				      // Read free Dirichlet multipliers
				      MQ=CreateMatrix(nstates,nstates);
				      if (!(err=ReadMatrix(f_in,idformat="//== Free Dirichlet multipliers %s",trailer,MQ)))
					if (sv=CreateMarkovStateVariable_Single(nstates,nobs,Prior,FreeDim,NonZeroIndex,MQ))
					  if (nlags_encoded > 0)
					    {
					      dw_FreeArray(sv->lag_index);
					      sv->nlags_encoded=nlags_encoded;
					      sv->nbasestates=nbasestates;
					      sv->lag_index=CreateLagIndex(sv->nbasestates,sv->nlags_encoded,sv->nstates);
					    }
				      FreeMatrix(MQ);
				    }
				  dw_FreeArray(NonZeroIndex);
				}
			      dw_FreeArray(FreeDim);
			    }
			}
		      FreeMatrix(Prior);

		      break;
		    }
		}
	    }
	}
    }
  if (err) ReadError(idformat,trailer,err);
  if (trailer) free(trailer);
  return sv;
}

int WriteMarkovSpecification_SV(FILE *f_out, TMarkovStateVariable *sv, char *idstring)
{
  int i, j;
  char *trailer, *idbuffer;

  if (idstring[0])
    {
      // 24 characters in "for state_variable ==//" plus null character
      trailer=(char*)malloc(24+strlen(idstring));
      sprintf(trailer,"for state_variable%s ==//",idstring);

      fprintf(f_out,"//****** Specification %s ******//\n\n",trailer);
    }
  else
    {
      trailer=(char*)malloc(5);
      strcpy(trailer,"==//");

      fprintf(f_out,"//== Number observations ==//\n%d\n\n",sv->nobs);
      fprintf(f_out,"//****** State Variable Specification ******//\n");
    }

  fprintf(f_out,"//== Number independent state variables %s\n%d\n",trailer,sv->n_state_variables);

  if (sv->n_state_variables > 1)
    {
      fprintf(f_out,"//******************************************//\n\n");
      for (j=10, i=1; sv->n_state_variables >= j; j*=10, i++);
      strcpy(idbuffer=(char*)malloc((j=(int)strlen(idstring))+i+3),idstring); 
      for (i=0; i < sv->n_state_variables; i++)
	{
	  sprintf(idbuffer+j,"[%d]",i+1);
	  if (!WriteMarkovSpecification_SV(f_out,sv->state_variable[i],idbuffer)) 
            {
	      free(idbuffer);
	      return 0;
	    }
	}
      free(idbuffer);
    }
  else
    {
      fprintf(f_out,"\n//== Number states %s\n%d\n\n",trailer,sv->nstates);

      if (sv->nlags_encoded > 0)
	{
	  fprintf(f_out,"//== Number of lags encoded %s\n%d\n\n",trailer,sv->nlags_encoded);
	  fprintf(f_out,"//== Number of base states %s\n%d\n\n",trailer,sv->nbasestates);
	}

      fprintf(f_out,"//== Prior %s\n",trailer);
      dw_PrintMatrix(f_out,sv->Prior,"%22.14le ");
      fprintf(f_out,"\n");

      fprintf(f_out,"//== Number free Dirichlet variables %s\n%d\n\n",trailer,dw_DimA(sv->FreeDim));

      fprintf(f_out,"//== Free Dirichlet dimensions %s\n",trailer);
      dw_PrintArray(f_out,sv->FreeDim,"%d ");
      fprintf(f_out,"\n");

      fprintf(f_out,"//== Free Dirichlet index %s\n",trailer);
      dw_PrintArray(f_out,sv->NonZeroIndex,"%d ");

      fprintf(f_out,"//== Free Dirichlet multipliers %s\n",trailer);
      dw_PrintMatrix(f_out,sv->MQ,"%22.14le ");
      fprintf(f_out,"//******************************************//\n\n");
    }

  free(trailer);
  return 1;
}

/*
   Assumes:
    f:  valid file pointer
    sv:  pointer to valid TMarkovStateVariable structure
    idstring:  pointer to null terminated string

   Returns:
    One upon success and zero upon failure.

   Results:
    Reads transition matrices from file into sv.  All transition matrices and the
    associated free parameters are set.  

   Notes:
    The null terminated string idstring is of the form "" or "[i][j]...[k]".  
    Usually, this function will be called with idstring equal to "".  The routine 
    is called recursively with the more general form for idstring.  For each path 
    of Markov state variables of the form

        sv->state_variable[i]->state_variable[j]-> ... ->state_variable[k]

    which corresponds to a single Markov state variable, there must be an entry 
    in the file of the form:

    //== <header>Transition matrix[i][j]...[k] ==//
    x x ... x
    x x ... x
    . .     .
    x x ... x

    If sv itself is a single Markov state variable, then the format can be 

    //== <header>Transition matrix[1] ==//
    x x ... x
    x x ... x
    . .     .
    x x ... x

    or

    //== <header>Transition matrix ==//
    x x ... x
    x x ... x
    . .     .
    x x ... x

    Here the term <header> is replaced with the null terminated string header.
    Note that the spacing is important.
*/
int ReadTransitionMatrices_SV(FILE *f_in, TMarkovStateVariable* sv, char *header, char *idstring)
{
  int i, j, err;
  char *format, *idbuffer;
  PRECISION sum;

  if (sv->n_state_variables > 1)
    {
      for (j=10, i=1; sv->n_state_variables >= j; j*=10, i++);
      strcpy(idbuffer=(char*)malloc((j=(int)strlen(idstring))+i+3),idstring); 
      for (i=sv->n_state_variables-1; i >= 0; i--)
	{
	  sprintf(idbuffer+j,"[%d]",i+1);
	  if (!ReadTransitionMatrices_SV(f_in,sv->state_variable[i],header,idbuffer)) 
            {
	      free(idbuffer);
	      return sv->valid_transition_matrix=0;
	    }
	}
      free(idbuffer);
      MatrixTensor(sv->Q,sv->QA);
      return sv->valid_transition_matrix=1;
    }
  else
    {
      // Read transition matrix
      if (!header) header="";
      format="//== %sTransition matrix%s ==//";
      sprintf(idbuffer=(char*)malloc(strlen(header) + strlen(format) + strlen(idstring) - 3),format,header,idstring);
      if (err=ReadMatrix(f_in,idbuffer,(char*)NULL,sv->Q))
	if (!idstring[0])
	  {
	    free(idbuffer);
	    idstring="[1]";
	    sprintf(idbuffer=(char*)malloc(strlen(header) + strlen(format) + strlen(idstring) - 3),format,header,idstring);
	    err=ReadMatrix(f_in,idbuffer,(char*)NULL,sv->Q);
	  }
      free(idbuffer);
      if (!err)
	{      
	  // Scale the columns of Q - loose requirement on sumation to one
	  for (j=sv->nstates-1; j >= 0; j--)
	    {
	      for (sum=0.0, i=sv->nstates-1; i >= 0; i--) 
		if (ElementM(sv->Q,i,j) < 0.0)
		  {
		    dw_UserError("Transition matrix can not have negative elements.");
		    return sv->valid_transition_matrix=0;
		  }
		else
		  sum+=ElementM(sv->Q,i,j);
	      if (fabs(sum-1.0) > 1.0e-4)
		{
		  dw_UserError("Transition matrix columns must sum to one.");
		  return sv->valid_transition_matrix=0;
		}
	      for (sum=1.0/sum, i=sv->nstates-1; i >= 0; i--) 
		ElementM(sv->Q,i,j)*=sum;
	    }

	  // Update
	  if (!Update_B_from_Q_SV(sv))
	    {
	      dw_UserError("Transition matrices do not satisfy restrictions");
	      return sv->valid_transition_matrix=0;
	    }

	  return sv->valid_transition_matrix=1;
	}
      return sv->valid_transition_matrix=0;
    }
}

/*
   Assumes:
    f:  valid file pointer
    sv:  pointer to valid TMarkovStateVariable structure
    idstring:  pointer to null terminated string

   Returns:
    One upon success and zero upon failure.

   Results:
    Reads transition matrices from file into sv.  All transition matrices and the
    associated free parameters are set.  

   Notes:
    The null terminated string idstring is of the form "" or "[i][j]...[k]".  
    Usually, this function will be called with idstring equal to "".  The routine 
    is called recursively with the more general form for idstring.  For each path 
    of Markov state variables of the form

        sv->state_variable[i]->state_variable[j]-> ... ->state_variable[k]

    which corresponds to a single Markov state variable, there must be an entry 
    in the file of the form:

    //== <header>Transition matrix[i][j]...[k] ==//
    x x ... x
    x x ... x
    . .     .
    x x ... x

    If sv itself is a single Markov state variable, then the format can be 

    //== <header>Transition matrix[1] ==//
    x x ... x
    x x ... x
    . .     .
    x x ... x

    or

    //== <header>Transition matrix ==//
    x x ... x
    x x ... x
    . .     .
    x x ... x

    Here the term <header> is replaced with the null terminated string header.
    Note that the spacing is important.
*/
int ReadBaseTransitionMatrices_SV(FILE *f_in, TMarkovStateVariable* sv, char *header, char *idstring)
{
  int i, j, err;
  char *format, *idbuffer;
  PRECISION sum;
  TMatrix Q;

  if (sv->n_state_variables > 1)
    {
      for (j=10, i=1; sv->n_state_variables >= j; j*=10, i++);
      strcpy(idbuffer=(char*)malloc((j=(int)strlen(idstring))+i+3),idstring); 
      for (i=sv->n_state_variables-1; i >= 0; i--)
	{
	  sprintf(idbuffer+j,"[%d]",i+1);
	  if (!ReadBaseTransitionMatrices_SV(f_in,sv->state_variable[i],header,idbuffer)) 
            {
	      free(idbuffer);
	      return 0;
	    }
	}
      free(idbuffer);
      MatrixTensor(sv->Q,sv->QA);
      return 1;
    }
  else
    {
      // Read transition matrix
      Q=CreateMatrix(sv->nbasestates,sv->nbasestates);
      if (!header) header="";
      format="//== %sBase transition matrix%s ==//";
      sprintf(idbuffer=(char*)malloc(strlen(header) + strlen(format) + strlen(idstring) - 3),format,header,idstring);
      if (err=ReadMatrix(f_in,idbuffer,(char*)NULL,Q))
	if (!idstring[0])
	  {
	    free(idbuffer);
	    idstring="[1]";
	    sprintf(idbuffer=(char*)malloc(strlen(header) + strlen(format) + strlen(idstring) - 3),format,header,idstring);
	    err=ReadMatrix(f_in,idbuffer,(char*)NULL,Q);
	  }
      free(idbuffer);
      if (!err)
	{      
	  // Scale the columns of Q - loose requirement on sumation to one
	  for (j=sv->nbasestates-1; j >= 0; j--)
	    {
	      for (sum=0.0, i=sv->nbasestates-1; i >= 0; i--) 
		if (ElementM(Q,i,j) < 0.0)
		  {
		    FreeMatrix(Q);
		    dw_UserError("Transition matrix can not have negative elements.");
		    return 0;
		  }
		else
		  sum+=ElementM(Q,i,j);
	      if (fabs(sum-1.0) > 1.0e-4)
		{
		  FreeMatrix(Q);
		  dw_UserError("Transition matrix columns must sum to one.");
		  return 0;
		}
	      for (sum=1.0/sum, i=sv->nbasestates-1; i >= 0; i--) 
		ElementM(Q,i,j)*=sum;
	    }

	  // Convert base transition matrix to full transition matrix.
          ConvertBaseTransitionMatrix(sv->Q,Q,sv->nlags_encoded);

	  // Update
	  if (!Update_B_from_Q_SV(sv))
	    {
	      dw_UserError("Transition matrices do not satisfy restrictions");
	      return 0;
	    }

	  return 1;
	}
      return 0;
    }
}


/*
   Assumes:
    f:  valid file pointer or null pointer
    sv:  pointer to valid TMarkovStateVariable structure
    header: pointer to null terminated string
    idstring:  pointer to null terminated string

   Returns:
    One is always returned.

   Results:
    Writes transition matrices from sv to a file.  See ReadTransitionMatrices()
    for the format.

   Notes:
    The null terminated string idstring is of the form  "" or "[i][j]...[k]".  
    Usually, this routine will be called with idstring equal to "".  The routine 
    is called recursively with the more general form for idstring. 
*/
int WriteTransitionMatrices_SV(FILE *f_out, TMarkovStateVariable* sv, char *header, char *idstring)
{
  int i, j;
  char *idbuffer;

  if (!header) header="";
  fprintf(f_out,"//== %sTransition matrix%s ==//\n",header,idstring);
  dw_PrintMatrix(f_out,sv->Q,"%22.14le ");

  if (sv->n_state_variables > 1)
    {
      for (j=10, i=1; sv->n_state_variables >= j; j*=10, i++);
      strcpy(idbuffer=(char*)malloc((j=(int)strlen(idstring))+i+3),idstring);
      for (i=0; i < sv->n_state_variables; i++)
	{
	  sprintf(idbuffer+j,"[%d]",i+1);
	  WriteTransitionMatrices_SV(f_out,sv->state_variable[i],header,idbuffer);
	}
      free(idbuffer);
    }
  
  return 1;
}

/*
   Assumes:
    f:  valid file pointer or null pointer
    sv:  pointer to valid TMarkovStateVariable structure
    header: pointer to null terminated string
    idstring:  pointer to null terminated string

   Returns:
    One is always returned.

   Results:
    Writes transition matrices from sv to a file.  See ReadTransitionMatrices()
    for the format.

   Notes:
    The null terminated string idstring is of the form  "" or "[i][j]...[k]".  
    Usually, this routine will be called with idstring equal to "".  The routine 
    is called recursively with the more general form for idstring. 
*/
int WriteBaseTransitionMatrices_SV(FILE *f_out, TMarkovStateVariable* sv, char *header, char *idstring)
{
  int i, j;
  char *idbuffer;
  TMatrix Q;

  if (!header) header="";
  fprintf(f_out,"//== %sBase transition matrix%s ==//\n",header,idstring);
  if (Q=GetBaseTransitionMatrix_SV((TMatrix)NULL,sv))
    {
      dw_PrintMatrix(f_out,Q,"%22.14le ");
      fprintf(f_out,"\n");
      FreeMatrix(Q);
    }
  else
    fprintf(f_out,"Error geting base transition matrix\n");

  if (sv->n_state_variables > 1)
    {
      for (j=10, i=1; sv->n_state_variables >= j; j*=10, i++);
      strcpy(idbuffer=(char*)malloc((j=(int)strlen(idstring))+i+3),idstring);
      for (i=0; i < sv->n_state_variables; i++)
	{
	  sprintf(idbuffer+j,"[%d]",i+1);
	  WriteBaseTransitionMatrices_SV(f_out,sv->state_variable[i],header,idbuffer);
	}
      free(idbuffer);
    }
  
  return 1;
}

/*
   Assumes:
    f:  valid file pointer or null pointer
    sv:  pointer to valid TMarkovStateVariable structure
    header: pointer to null terminated string
    idstring:  pointer to null terminated string

   Returns:
    One is always returned.

   Results:
    Writes transition matrices from sv to a file.  See ReadTransitionMatrices()
    for the format.

   Notes:
    The null terminated string idstring is of the form  "" or "[i][j]...[k]".  
    Usually, this routine will be called with idstring equal to "".  The routine 
    is called recursively with the more general form for idstring. 
*/
void WriteBaseTransitionMatricesFlat_Headers_SV(FILE *f_out, TMarkovStateVariable* sv, char *idstring)
{
  int i, j;
  char *idbuffer;

  if (sv->n_state_variables > 1)
    {
      for (j=10, i=1; sv->n_state_variables >= j; j*=10, i++);
      strcpy(idbuffer=(char*)malloc((j=(int)strlen(idstring))+i+3),idstring);
      for (i=0; i < sv->n_state_variables; i++)
	{
	  sprintf(idbuffer+j,"[%d]",i+1);
	  WriteBaseTransitionMatricesFlat_Headers_SV(f_out,sv->state_variable[i],idbuffer);
	}
      free(idbuffer);
    }
  else
    {
      for (j=0; j < sv->nbasestates; j++)
	for (i=0; i < sv->nbasestates; i++)
	  fprintf(f_out,"Q%s(%d,%d) ",idstring,i+1,j+1);
    }
}


/*
   Returns 1 upon success and 0 upon failure.
*/
int WriteBaseTransitionMatricesFlat_SV(FILE *f_out, TMarkovStateVariable *sv, char *fmt)
{
  int i, j;
  TMatrix Q;

  if (sv->n_state_variables > 1)
    {
      for (i=0; i < sv->n_state_variables; i++)
	if (!WriteBaseTransitionMatricesFlat_SV(f_out,sv->state_variable[i],fmt))
	  return 0;
    }
  else
    {
      if (!fmt) fmt="%lf ";
      if (Q=GetBaseTransitionMatrix_SV((TMatrix)NULL,sv))
	{
	  for (j=0; j < ColM(Q); j++)
	    for (i=0; i < RowM(Q); i++)
	      fprintf(f_out,fmt,ElementM(Q,i,j));
	  FreeMatrix(Q);
	}
      else
	return 0;
    }

  return 1;
}
#undef SWITCHIO_LINE_ID_NOT_FOUND 
#undef SWITCHIO_ERROR_READING_DATA
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*
   Assumes:
    f:  valid file pointer or null pointer
    filename:  pointer to null terminated string or null

   Returns:
    Valid pointer to a TMarkovStateVariable structure upon success and null 
    pointer upon failure.

   Results:
    Reads Markov specification from file and creates TMarkovStateVariable 
    structure. 

   Notes:
    One of f or filename should be non-null.
*/
TMarkovStateVariable* ReadMarkovSpecification(FILE *f, char *filename)
{
  TMarkovStateVariable *sv;
  FILE *f_in=f ? f : fopen(filename,"rt");
  sv=(f_in) ? ReadMarkovSpecification_SV(f_in,"",-1) : (TMarkovStateVariable*)NULL;
  if (!f && f_in) fclose(f_in);
  return sv;
}

/*
   Assumes:
    f:  valid file pointer or null pointer
    filename:  pointer to null terminated string or null
    model: pointer to valid TStateModel structure.

   Returns:
    one upon upon success and zero upon failure.

   Results:
    Writes Markov specification to file. 

   Notes:
    One of f or filename should be non-null.
*/
int WriteMarkovSpecification(FILE *f, char *filename, TStateModel *model)
{
  FILE *f_out=f ? f : fopen(filename,"at");
  int rtrn=(f_out) ? WriteMarkovSpecification_SV(f_out,model->sv,"") : 0;
  if (!f && f_out) fclose(f_out);
  return rtrn;
}

/*
   Assumes:
    f:  valid file pointer or null pointer
    filename:  pointer to null terminated string or null
    model: pointer to valid TStateModel structure.

   Returns:
    One upon success and zero upon failure.

   Results:
    Attempts to read transition matrix into model.

   Notes:
    One of f or filename should be non-null.  Calls TransitionMatricesChanged().
*/
int ReadTransitionMatrices(FILE *f, char *filename, char *header, TStateModel *model)
{
  FILE *f_in=f ? f : fopen(filename,"rt");
  int rtrn=(f_in) ? ReadTransitionMatrices_SV(f_in,model->sv,header,"") : 0;
  TransitionMatricesChanged(model);
  if (!f && f_in) fclose(f_in);
  return rtrn;
}

/*
   Assumes:
    f:  valid file pointer or null pointer
    filename:  pointer to null terminated string or null
    model: pointer to valid TStateModel structure.

   Returns:
    One upon success and zero upon failure.

   Results:
    Writes transition matrix.

   Notes:
    One of f or filename should be non-null. 
*/
int WriteTransitionMatrices(FILE *f, char *filename, char *header, TStateModel *model)
{
  FILE *f_out=f ? f : dw_CreateTextFile(filename);
  int rtrn=(f_out) ? WriteTransitionMatrices_SV(f_out,model->sv,header,"") : 0;
  if (!f && f_out) fclose(f_out);
  return rtrn;
}

/*
   Assumes:
    f:  valid file pointer or null pointer
    filename:  pointer to null terminated string or null
    model: pointer to valid TStateModel structure.

   Returns:
    One upon success and zero upon failure.

   Results:
    Attempts to read transition matrix into model.

   Notes:
    One of f or filename should be non-null.  Calls TransitionMatricesChanged().
*/
int ReadBaseTransitionMatrices(FILE *f, char *filename, char *header, TStateModel *model)
{
  FILE *f_in=f ? f : fopen(filename,"rt");
  int rtrn=(f_in) ? ReadBaseTransitionMatrices_SV(f_in,model->sv,header,"") : 0;
  TransitionMatricesChanged(model);
  if (!f && f_in) fclose(f_in);
  return rtrn;
}

/*
   Assumes:
    f:  valid file pointer or null pointer
    filename:  pointer to null terminated string or null
    model: pointer to valid TStateModel structure.

   Returns:
    One upon success and zero upon failure.

   Results:
    Writes base transition matrices.

   Notes:
    One of f or filename should be non-null. 
*/
int WriteBaseTransitionMatrices(FILE *f, char *filename, char *header, TStateModel *model)
{
  FILE *f_out=f ? f : dw_CreateTextFile(filename);
  int rtrn=(f_out) ? WriteBaseTransitionMatrices_SV(f_out,model->sv,header,"") : 0;
  if (!f && f_out) fclose(f_out);
  return rtrn;
}

int WriteBaseTransitionMatricesFlat(FILE *f, TStateModel *model, char *fmt)
{
  return f ? WriteBaseTransitionMatricesFlat_SV(f,model->sv,fmt) : 0;
}


int ReadStates(FILE *f, char *filename, char *header, TStateModel *model)
{
  FILE *f_in=f ? f : dw_OpenTextFile(filename);
  char *format="//== %sStates ==//";
  int err, i;

  if (err=ReadIntArray(f_in,format,header,model->sv->S))
    ReadError(format,header,err);
  else
    {
      // Check states and propagate
      for (i=model->sv->nstates; i >= 0; i--)
	if ((model->sv->S[i] < 0) || (model->sv->S[i] >= model->sv->nstates))
	  {
	    for ( ; i >= 0; i--) model->sv->S[i]=0;
	    dw_UserError("ReadStates(): Invalid state value.");
	    err=1;
	    break;
	  }
      PropagateStates_SV(model->sv);
    }

  if (!f) fclose(f_in);

  return err ? 0 : 1;
}

int WriteStates(FILE *f, char *filename, char *header, TStateModel *model)
{
  FILE *f_out=f ? f : dw_CreateTextFile(filename);

  fprintf(f_out,"//== %sStates ==//\n",header); 
  dw_PrintArray(f_out,model->sv->S,"%d ");
  fprintf(f_out,"\n");

  if (!f) fclose(f_out);
  return 1;
}

/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
/*
//== Flat Independent Markov States and Simple Restrictions ==//

//-----------------------------------------------------------------------------//
//-- Read by CreateMarkovStateVariable_File() only if the passed number of   --//
//-- observations is less than or equal to zero.  Can be omitted if the      --//
//-- passed number of observations is positive.                              --//
//-----------------------------------------------------------------------------//
//== Number Observations ==//
nobs

//== Number Independent State Variables ==//
n_state_variables

//-----------------------------------------------------------------------------//
//-- state_variable[i]                                                       --//
//-----------------------------------------------------------------------------//

//== Number of states for state_variable[i] ==//
n_states

//-----------------------------------------------------------------------------//
//-- Each column contains the parameters for a Dirichlet prior on the        --//
//-- corresponding column of the transition matrix.  Each element must be    --//
//-- positive.  For each column, the relative size of the prior elements     --//
//-- determine the relative size of the elements of the transition matrix    --//
//-- and overall larger sizes implies a tighter prior.                       --//
//-----------------------------------------------------------------------------//
//== Transition matrix prior for state_variable[i]. (n_states x n_states) ==//
prior

//== Free Dirichet dimensions for state_variable[i]  ==//
free[0] ... free[n_states - 1]

//-----------------------------------------------------------------------------//
//-- The jth restriction matrix is n_states x free[j].  Each row of the      --//
//-- restriction matrix has exactly one non-zero entry and the sum of each   --//
//-- column of the restriction matrix must be one.  If entry (i,k) of the    --//
//-- jth restriction matrix is non-zero, then entry (i,j) in the transition  --//
//-- matrix is controlled by the kth element of jth free Dirichlet random    --//
//-- variable                                                                --//
//-----------------------------------------------------------------------------//
//== Column restrictions for state_variable[i] ==//
restriction[0]
     .
     .
     .
restriction[n_states - 1]

//-----------------------------------------------------------------------------//
//-- Allows an optional initialization for the transition matrix, otherwise  --//
//-- the initial value is set to the mean of the prior                       --//
//-----------------------------------------------------------------------------//
//== Initial transition matrix for state_variable[i]. (n_states x n_states)  ==//
Q

//-----------------------------------------------------------------------------//
//-- Allows for lagged values of the state variable to be encoded.  If this  --//
//-- identifier is missing, then the value of nlags_encoded is set to zero.  --//
//-----------------------------------------------------------------------------//
//== Number of lags encoded for state_variable[i] ==//
nlags_encoded

*/
TMarkovStateVariable* CreateMarkovStateVariable_File(FILE *f, char *filename, int nobs)
{
  FILE *f_in;
  char *id, id_buffer[256];
  int nlags, n_state_variables, nstates, i, j;
  int* dims;
  TMatrix prior;
  TMatrix* restrictions;
  TMarkovStateVariable **sv, *rtrn=(TMarkovStateVariable*)NULL, *tmp;

  // Open file if necessary
  if (!f)
    f_in=dw_OpenTextFile(filename);
  else
    f_in=f;

  // Check for Flat Independent Markov States and Simple Restrictions
  if (dw_SetFilePosition(f_in,"//== Flat Independent Markov States and Simple Restrictions ==//"))
    {
      if (nobs <= 0)
	{
	  id="//== Number Observations ==//";
	  if (!dw_SetFilePosition(f_in,id))
	    {
	      fprintf(stderr,"Line identifier ""%s"" not found.\n",id);
	      exit(0);
	    }
	  fscanf(f_in," %d ",&nobs);
	  if (nobs <= 0)
	    {
	      fprintf(stderr,"Number Observations must be positive\n");
	      exit(0);
	    }
	}

      id="//== Number Independent State Variables ==//";
      if (!dw_SetFilePosition(f_in,id))
        {
          fprintf(stderr,"Line identifier ""%s"" not found.\n",id);
          exit(0);
        }
      fscanf(f_in," %d ",&n_state_variables);
      if (n_state_variables <= 0)
	{
	  fprintf(stderr,"Number Independent State Variables must be positive\n");
	  exit(0);
	}

      sv=(TMarkovStateVariable**)dw_CreateArray_pointer(n_state_variables,(void (*)(void*))FreeMarkovStateVariable);
      for (i=0; i < n_state_variables; i++)
        {
	  sprintf(id_buffer,"//== Number of states for state_variable[%d] ==//",i+1);
	  if (!dw_SetFilePosition(f_in,id_buffer))
	    {
	      fprintf(stderr,"Line identifier ""%s"" not found.\n",id_buffer);
	      exit(0);
	    }
	  fscanf(f_in," %d ",&nstates);
	  if (nstates <= 0)
	    {
	      fprintf(stderr,"Number of states for state_variable[%d] must be positive\n",i+1);
	      exit(0);
	    }

	  sprintf(id_buffer,"//== Transition matrix prior for state_variable[%d]. (n_states x n_states) ==//",i+1);
	  if (!dw_SetFilePosition(f_in,id_buffer))
	    {
	      fprintf(stderr,"Line identifier ""%s"" not found.\n",id_buffer);
	      exit(0);
	    }
	  dw_ReadMatrix(f_in,prior=CreateMatrix(nstates,nstates));

	  sprintf(id_buffer,"//== Free Dirichet dimensions for state_variable[%d]  ==//",i+1);
	  if (!dw_SetFilePosition(f_in,id_buffer))
	    sv[i]=CreateMarkovStateVariable_NoRestrictions(nstates,nobs,prior);
	  else
	    {
	      dw_ReadArray(f_in,dims=dw_CreateArray_int(nstates));

	      sprintf(id_buffer,"//== Column restrictions for state_variable[%d] ==//",i+1);
	      if (!dw_SetFilePosition(f_in,id_buffer))
		{
		  fprintf(stderr,"Line identifier ""%s"" not found.\n",id_buffer);
		  exit(0);
		}
	      restrictions=dw_CreateArray_matrix(nstates);
	      for (j=0; j < nstates; j++)
		if (dims[j] > 0)
		  dw_ReadMatrix(f_in,restrictions[j]=CreateMatrix(nstates,dims[j]));
		else
		  {
		    fprintf(stderr,"Free Dirichet dimensions for column %d of state_variable[%d] must be positive\n",j+1,i+1);
		    exit(0);
		  }

	      sv[i]=CreateMarkovStateVariable_SimpleRestrictions(nstates,nobs,prior,restrictions);

	      dw_FreeArray(restrictions);
	      dw_FreeArray(dims);

              sprintf(id_buffer,"//== Number of lags encoded for state_variable[%d] ==//",i+1);
	      if (dw_SetFilePosition(f_in,id_buffer))
		{
		  fscanf(f_in," %d ",&nlags);
		  if (nlags > 0)
		    {
		      tmp=CreateMarkovStateVariable_Lags(nlags,sv[i]);
		      FreeMarkovStateVariable(sv[i]);
		      sv[i]=tmp;
		    }
		}
	    }

	  FreeMatrix(prior);
	}

      if (n_state_variables > 1)
	rtrn=CreateMarkovStateVariable_Multiple(nobs,n_state_variables,sv);
      else
	rtrn=sv[0];
    }

  // Close file if necessary
  if (!f) fclose(f_in);

  return rtrn;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


/*******************************************************************************/
/************************* Test Input/Output Routines **************************/
/*******************************************************************************
int main(void)
{
  char *filename="../switch/MarkovStateVariable.dat";
  char *outputfilename="../switch/tmp.dat";
  FILE *f_out;

  TMarkovStateVariable *sv;

  sv=CreateMarkovStateVariable_File((FILE*)NULL,filename,-1);

  f_out=dw_CreateTextFile(outputfilename);
  WriteMarkovSpecification_SV(f_out,sv,"");
  WriteTransitionMatrices_SV(f_out,sv,"Initial: ","");

  //DrawStatesFromTransitionMatrix_SV(sv);
  DrawTransitionMatrixFromPrior_SV(sv);
  WriteTransitionMatrices_SV(f_out,sv,"Draw: ","");
  fclose(f_out);

  FreeMarkovStateVariable(sv);

  return 0;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
