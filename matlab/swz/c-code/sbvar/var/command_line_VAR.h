#ifndef __command_line_VAR
#define __command_line_VAR

#include "switch.h"
#include <stdio.h>

char* CreateFilenameFromTag(char *fmt, char *tag, char *dir);
char* CreatePath(char *path);
void PrintHelpMessages(FILE *f, char **include, char **additional);

typedef struct
{
  char *in_directory;              // -di
  char *in_tag;                    // -ft
  char *specification_filename;    // -fs
  char *parameters_filename;       // -fp
  char *parameters_header;         // -ph

  char *specification_filename_actual;
  char *parameters_filename_actual;
  char *parameters_header_actual;

  char *out_directory;             // -do
  char *out_tag;                   // -fto (default from -ft)
  char *out_header;                // -pho (default from -ph)
} TVARCommandLine;

TVARCommandLine* Create_VARCommandLine(void);
void Free_VARCommandLine(TVARCommandLine *cmd);
TVARCommandLine* Base_VARCommandLine(int nargs, char **args, TVARCommandLine *cmd);

void EstimateFinal_VARCommandLine_Help(FILE *f);
TStateModel* CreateTStateModelFromEstimateFinal(int nargs, char **args, TVARCommandLine **p_cmd);

TStateModel* CreateTStateModelForEstimate(int nargs, char **args, TVARCommandLine **p_cmd);

#endif
