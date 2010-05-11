
#ifndef __PARSE_COMMAND_LINE__
#define __PARSE_COMMAND_LINE__

#ifdef __cplusplus
extern "C"
{
#endif

int dw_FindArgument(int nargs, char **args, char opt);
int dw_ParseInteger(int nargs, char **args, char opt, int def);
double dw_ParseFloating(int nargs, char **args, char opt, double def);
char* dw_ParseString(int nargs, char **args, char opt, char *def);

int dw_FindArgument_String(int nargs, char **args, char *opt);
int dw_ParseInteger_String(int nargs, char **args, char *opt, int def);
double dw_ParseFloating_String(int nargs, char **args, char *opt, double def);
char* dw_ParseString_String(int nargs, char **args, char *opt, char *def);

#ifdef __cplusplus
}
#endif

#endif
