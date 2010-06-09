
#include "dw_error.h"
#include <stdlib.h>
#include <string.h>

#include "modify_for_mex.h"

#define ERROR_MESSAGE_BUFFER_LENGTH 256
static int ERROR=NO_ERR;
static char ERROR_MESSAGE[ERROR_MESSAGE_BUFFER_LENGTH]="";
static int TerminalErrors=ALL_ERRORS;
static int VerboseErrors=ALL_ERRORS;
static FILE* f_err=(FILE*)NULL;

/*
   Returns the value of the current error flag.
*/
int dw_GetError(void)
{
  return ERROR;
}

/*
   Returns pointer to current error message.  This buffer should not be modified
   or freed.
*/
char* dw_GetErrorMessage(void)
{
  return ERROR_MESSAGE;
}

/*
   Clears the error flag and returns the value of the previous flag.  This is the
   most efficient way to clear the error flag and message.
*/
int dw_ClearError(void)
{
  int rtrn=ERROR;
  ERROR=NO_ERR;
  ERROR_MESSAGE[0]='\0';
  return rtrn;
}


/*
   Sets the error flag to err, the error message to msg and returns the value of
   the previous flag.  If the null terminated string msg is longer than 255
   characters, only the first 255 characters are used.  If msg is null, then a
   predefined error message is used.
*/
int dw_SetError(int err, char *msg)
{
  int rtrn=ERROR;
  if (msg)
    switch (ERROR=err)
      {
      case MEM_ERR:
      case FILE_ERR:
      case PARSE_ERR:
      case FLOAT_ERR:
      case NULL_ERR:
      case ARG_ERR:
      case ITERATION_ERR:
      case NOT_IMPLEMENTED_ERR:
      case SIZE_ERR:
      case SING_ERR:
      case POSDEF_ERR:
      case BLAS_LAPACK_ERR:
      case USER_ERR:
    strncpy(ERROR_MESSAGE,msg,ERROR_MESSAGE_BUFFER_LENGTH-1);
    ERROR_MESSAGE[ERROR_MESSAGE_BUFFER_LENGTH-1]='\0';
    break;
      case NO_ERR:
    ERROR_MESSAGE[0]='\0';
    break;
      default:
    ERROR=UNKNOWN_ERR;
    strcpy(ERROR_MESSAGE,"Unknown error.");
    break;
      }
  else
    switch (ERROR=err)
      {
      case MEM_ERR:
    strcpy(ERROR_MESSAGE,"Out of memory.");
    break;
      case FILE_ERR:
    strcpy(ERROR_MESSAGE,"File operation error.");
    break;
      case PARSE_ERR:
    strcpy(ERROR_MESSAGE,"Error parsing data.");
    break;
      case FLOAT_ERR:
    strcpy(ERROR_MESSAGE,"Floating point error.");
    break;
      case NULL_ERR:
    strcpy(ERROR_MESSAGE,"Unexpected null pointer encountered.");
    break;
      case ARG_ERR:
    strcpy(ERROR_MESSAGE,"Argument error.");
    break;
      case ITERATION_ERR:
    strcpy(ERROR_MESSAGE,"Maximum iteration limit exceeded.");
    break;
      case NOT_IMPLEMENTED_ERR:
    strcpy(ERROR_MESSAGE,"Feature not yet implemented.");
    break;
      case SIZE_ERR:
    strcpy(ERROR_MESSAGE,"Matrices/vectors not conformable.");
    break;
      case SING_ERR:
    strcpy(ERROR_MESSAGE,"Singular matrix.");
    break;
      case POSDEF_ERR:
    strcpy(ERROR_MESSAGE,"Matrix not positive definite.");
    break;
      case BLAS_LAPACK_ERR:
    strcpy(ERROR_MESSAGE,"Blas/Lapack error.");
    break;
      case USER_ERR:
        strcpy(ERROR_MESSAGE,"Undocumented error.");
    break;
      case NO_ERR:
    ERROR_MESSAGE[0]='\0';
    break;
      default:
    ERROR=UNKNOWN_ERR;
    strcpy(ERROR_MESSAGE,"Unknown error.");
    break;
      }
  if (VerboseErrors & ERROR) fprintf(f_err ? f_err : stderr,"%s\n",ERROR_MESSAGE);
  if (TerminalErrors & ERROR) exit(ERROR);
  return rtrn;
}

/*
   Sets the error flag and to err, sets the error message to the predefined error
   message, and returns the value of the previous error flag.
*/
int dw_Error(int err)
{
  return dw_SetError(err,(char*)NULL);
}

/*
   Sets the error flag and to USER_ERR, sets the error message to msg, and
   returns the value of the previous error flag.
*/
int dw_UserError(char *msg)
{
  return dw_SetError(USER_ERR,msg);
}

/*
   Sets errors which terminate program.  The integer err should be a combination
   of the error flags defined in dw_error.h.
*/
int dw_SetTerminalErrors(int err)
{
  int rtrn=TerminalErrors;
  TerminalErrors=err & ALL_ERRORS;
  return rtrn;
}

/*
  Returns the current terminal errors.
*/
int dw_GetTerminalErrors(void)
{
  return TerminalErrors;
}

/*
   Sets errors which causes program to print a error message to f_err.  The
   integer err should be a combination of the error flags defined in dw_error.h.
*/
int dw_SetVerboseErrors(int err)
{
  int rtrn=VerboseErrors;
  VerboseErrors=err & ALL_ERRORS;
  return rtrn;
}

/*
  Returns the current verbose errors.
*/
int dw_GetVerboseErrors(void)
{
  return VerboseErrors;
}

/*
   Sets the file to which errors messages will be sent.  The  file pointer f
   must either be the null pointer or a valid printer to an open file.  Passing a
   null pointer has the same effect as redirecting output to stderr.  To suppress
   output of error messages, call dw_SetVerboseErrors().  Returns a pointer to
   current error message file.  When redirecting the error message, if is
   critical that the error message file pointer point to an open file and that
   this file not be closed as long as it is the current error message file.
*/
FILE* dw_SetErrorMessageFile(FILE *f)
{
  FILE *rtrn=f_err;
  f_err=f ? f : stderr;
  return rtrn;
}


