#ifndef __command_line_VAR
#define __command_line_VAR

#include "switch.h"
#include <stdio.h>

char* CreateFilenameFromTag(char *fmt, char *tag, char *dir);
char* CreatePath(char *path);
void PrintHelpMessages(FILE *f, char **include, char **additional);

typedef struct
{
  char *in_directory;               /*   -di   ansi-c*/
  char *in_tag;                     /*   -ft   ansi-c*/
  char *specification_filename;     /*   -fs   ansi-c*/
  char *parameters_filename;        /*   -fp   ansi-c*/
  char *parameters_header;          /*   -ph   ansi-c*/

  char *specification_filename_actual;
  char *parameters_filename_actual;
  char *parameters_header_actual;

  char *out_directory;              /*   -do   ansi-c*/
  char *out_tag;                    /*   -fto (default from -ft)   ansi-c*/
  char *out_header;                 /*   -pho (default from -ph)   ansi-c*/
} TVARCommandLine;

TVARCommandLine* Create_VARCommandLine(void);
void Free_VARCommandLine(TVARCommandLine *cmd);
TVARCommandLine* Base_VARCommandLine(int nargs, char **args, TVARCommandLine *cmd);

void EstimateFinal_VARCommandLine_Help(FILE *f);
TStateModel* CreateTStateModelFromEstimateFinal(int nargs, char **args, TVARCommandLine **p_cmd);

TStateModel* CreateTStateModelForEstimate(int nargs, char **args, TVARCommandLine **p_cmd);

#endif
