
#ifndef __ERROR_HANDLING__
#define __ERROR_HANDLING__

#include <stdio.h>

#define NO_ERR                0x00000000
#define ALL_ERRORS            0x000F03FF

/*  //=== General Errors ===   ansi-c*/
#define MEM_ERR               0x00000001
#define FILE_ERR              0x00000002
#define PARSE_ERR             0x00000004
#define FLOAT_ERR             0x00000008
#define NULL_ERR              0x00000010
#define ARG_ERR               0x00000020
#define ITERATION_ERR         0x00000040
#define USER_ERR              0x00000080
#define NOT_IMPLEMENTED_ERR   0x00000100
#define UNKNOWN_ERR           0x00000200


/*  //=== Matrix Errors ===   ansi-c*/
#define SIZE_ERR              0x00010000
#define SING_ERR              0x00020000
#define POSDEF_ERR            0x00040000
#define BLAS_LAPACK_ERR       0x00080000

/*  //=== Error Routines ===   ansi-c*/
int dw_GetError(void);
char* dw_GetErrorMessage(void);
int dw_ClearError(void);
int dw_SetError(int err, char *msg);
int dw_Error(int err);
int dw_UserError(char *msg);
int dw_SetVerboseErrors(int errors);
int dw_GetVerboseErrors(void);
int dw_SetTerminalErrors(int errors);
int dw_GetTerminalErrors(void);
FILE* dw_SetMessageFile(FILE *f);

#endif
