

#ifndef __DW_ASCII_ROUTINES__
#define __DW_ASCII_ROUTINES__

#include <stdio.h>

/*  // Flag codes.  See ParseDelimitedString() for explanation.   ansi-c*/
#define REMOVE_EMPTY_FIELDS           0x00000001
#define ALLOW_QUOTED_TEXT             0x00000002
#define STRIP_LEADING_WHITESPACE      0x00000004
#define STRIP_TRAILING_WHITESPACE     0x00000008
#define STRIP_WHITESPACE              0x0000000c

FILE *dw_OpenTextFile(char *filename);
FILE *dw_CreateTextFile(char *filename);
FILE *dw_AppendTextFile(char *filename);

int dw_SetFilePosition(FILE *f, char *id);
int dw_SetFilePositionNoRewind(FILE *f, char *id);
int dw_SetFilePositionBySection(FILE *f, int n, ...);

char* dw_ReadLine(FILE *f, char *buffer, int *n);
char** dw_ParseDelimitedString(char *buffer, char delimiter, int flag);
char** dw_ReadDelimitedLine(FILE *f, char delimiter, int flag);
char*** dw_ReadDelimitedFile(FILE *f, char* filename, char delimiter, int flag);
int dw_PrintDelimitedArray(FILE *f, void* array, char delimiter);

/*  //int dw_ReadDelimitedField(FILE *f, char **buffer, int *n);   ansi-c*/
int dw_ReadDelimitedField(FILE *f, int delimiter, int terminal, int flag, char **buffer, int *n);

char* dw_DuplicateString(char *buffer);

#endif
