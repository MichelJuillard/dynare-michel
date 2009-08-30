
#include "dw_parse_cmd.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define ARGUMENT_ID '-'

/*
  A floating point number is of the form

  [white space][+/-]digits[.[digits]][E/e[+/-]digits]white space/null character

  or

  [white space][+/-].digits[E/e[+/-]digits]white space/null character

  where characters in square brackets are optional.

  Returns one if valid floating point number and zero otherwise.
*/
int dw_IsFloat(char *buffer)
{
 int i=0;

 if (!buffer) return 0;

 /* Strip leading white space */
 while (isspace(buffer[i])) i++;

 /* Mantissa OK? */
 if ((buffer[i] == '+') || (buffer[i] == '-')) i++;
 if (isdigit(buffer[i]))
   {
    while (isdigit(buffer[++i]));
    if ((buffer[i] == '.'))
     while (isdigit(buffer[++i]));
   }
  else
   if ((buffer[i] == '.'))
     if (isdigit(buffer[++i]))
       while (isdigit(buffer[++i]));
      else
       return 0;
    else
     return 0;

 /* Is exponent OK? */
 if ((buffer[i] == 'e') || (buffer[i] == 'E'))
  {
   if ((buffer[++i] == '+') || (buffer[i] == '-')) i++;
   if (isdigit(buffer[i]))
     while (isdigit(buffer[++i]));
    else
     return 0;
  }

 /* Is end of string or trailing white space */
 if (buffer[i] && !isspace(buffer[i])) return 0;

 return 1;
}

/*
  Integers are of the form

   [white space][+/-]digits[.]white space/null character

  where characters in square brackets are optional.

  Returns one if valid integer and zero otherwise.
*/
int dw_IsInteger(char *buffer)
{
 int i=0;

 if (!buffer) return 0;

 /* Strip leading white space */
 while (isspace(buffer[i])) i++;

 /* Leading sign */
 if ((buffer[i] == '+') || (buffer[i] == '-')) i++;

 /* At least one digits possibly followed by decimal point */
 if (isdigit(buffer[i]))
   {
    while (isdigit(buffer[++i]));
    if ((buffer[i] == '.')) i++;
   }
  else
   return 0;

 /* Is end of string or trailing white space */
 if (buffer[i] && !isspace(buffer[i])) return 0;

 return 1;
}

/* 
   Searches args for a leading ARGUMENT_ID followed by the character opt.  Returns
   the index if found and -1 otherwise.
*/
int dw_FindArgument(int nargs, char **args, char opt)
{
  int i;
  for (i=nargs-1; i >= 0; i--)
    if ((args[i][0] == ARGUMENT_ID) && (args[i][1] == opt)) break;
  return i;
}

/*
   Searches for the last argument whose leading character is ARGUMENT_ID 
   followed by the character opt.  If such an argument is not found, then the 
   integer def is returned.  If such an argument is found then:

   Case 1:  The string length of the found argument is greater than 2.
     If the characters following the second form a valid integer, then this
     integer is returned.  Otherwise the integer def is returned.

   Case 2:  The string length of the found argument is equal to 2. 
     If there is an i+1 argument and its characters form a valid integer, then
     this integer is returned.  Otherwise the integer def is returned.  
*/
int dw_ParseInteger(int nargs, char **args, char opt, int def)
{
  int i=dw_FindArgument(nargs,args,opt);
  if (i != -1)
    if (dw_IsInteger(args[i]+2))
      return atoi(args[i]+2);
    else
      if ((i+1 < nargs) && dw_IsInteger(args[i+1])) return atoi(args[i+1]);
  return def;
}

/*
   Searches for the last argument whose leading character is ARGUMENT_ID
   followed by the character opt.  If such an argument is not found, then the 
   double def is returned.  If such an argument is found then:

   Case 1:  The string length of the found argument is greater than 2.
     If the characters following the second form a valid floating point number, 
     then this value is returned.  Otherwise def is returned.

   Case 2:  The string length of the found argument is equal to 2. 
     If there is an i+1 argument and its characters form a valid floating point
     number, then this value is returned.  Otherwise def is returned.  
*/
double dw_ParseFloating(int nargs, char **args, char opt, double def)
{
  int i=dw_FindArgument(nargs,args,opt);
  if (i != -1)
    if (dw_IsFloat(args[i]+2))
      return atof(args[i]+2);
    else
      if ((i+1 < nargs) && dw_IsFloat(args[i+1])) return atof(args[i+1]);
  return def;
}


/*
   Searches for the last argument whose leading character is ARGUMENT_ID
   followed by the character opt.  If such an argument is not found, then the 
   pointer def is returned.  If such an argument is found then:

   Case 1:  The string length of the found argument is greater than 2.
     A pointer to the found argument plus two is returned.

   Case 2:  The string length of the found argument is equal to 2. 
     If there is an i+1 argument then a pointer to this argument is returned.
     Otherwise the integer def is returned.  
*/
char* dw_ParseString(int nargs, char **args, char opt, char *def)
{
  int i=dw_FindArgument(nargs,args,opt);
  if (i != -1)
    if (args[i][2])
      return args[i]+2;
    else
      if (i+1 < nargs) return args[i+1];
  return def;
}

/* 
   Searches args for a leading ARGUMENT_ID followed by the string opt.  Returns
   the index if found and -1 otherwise.
*/
int dw_FindArgument_String(int nargs, char **args, char *opt)
{
  int i;
  for (i=nargs-1; i >= 0; i--)
    if ((args[i][0] == ARGUMENT_ID) && !strcmp(args[i]+1,opt)) break;
  return i;
}

/*
   Searches for the last argument whose leading character is a ARGUMENT_ID
   followed by the string opt.  If such an argument is not found, then the 
   integer def is returned.  If such an argument is found then:

   Case 1:  The string length of the found argument is greater than 1+strlen(opt).  
     If the characters following the second form a valid integer, then this 
     integer is returned.  Otherwise the integer def is returned.

   Case 2:  The string length of the found argument is equal to 1+strlen(opt). 
     If there is an i+1 argument and its characters form a valid integer, then
     this integer is returned.  Otherwise the integer def is returned.  
*/
int dw_ParseInteger_String(int nargs, char **args, char *opt, int def)
{
  int i=dw_FindArgument_String(nargs,args,opt);
  if ((i != -1) && (i+1 < nargs) && dw_IsInteger(args[i+1])) return atoi(args[i+1]);
  return def;
}

/*
   Searches for the last argument whose leading character is ARGUMENT_ID
   followed by the string opt.  If such an argument is not found, then the 
   double def is returned.  If such an argument is found then:

   Case 1:  The string length of the found argument is greater than 1+strlen(opt).
     If the characters following the second form a valid floating point number, 
     then this value is returned.  Otherwise def is returned.

   Case 2:  The string length of the found argument is equal to 1+strlen(opt).
     If there is an i+1 argument and its characters form a valid floating point
     number, then this value is returned.  Otherwise def is returned.  
*/
double dw_ParseFloating_String(int nargs, char **args, char *opt, double def)
{
  int i=dw_FindArgument_String(nargs,args,opt);
  if ((i != -1) && (i+1 < nargs) && dw_IsFloat(args[i+1])) return atof(args[i+1]);
  return def;
}


/*
   Searches for the last argument whose leading character is ARGUMENT_ID
   followed by the string opt.  If such an argument is not found, then the 
   pointer def is returned.  If such an argument is found then:

   Case 1:  The string length of the found argument is greater than 1+strlen(opt).
     A pointer to the found argument plus two is returned.

   Case 2:  The string length of the found argument is equal to 1+strlen(opt). 
     If there is an i+1 argument, then a pointer to this argument is returned.  
     Otherwise the string def is returned.  
*/
char* dw_ParseString_String(int nargs, char **args, char *opt, char *def)
{
  int i=dw_FindArgument_String(nargs,args,opt);
  if ((i != -1) && (i+1 < nargs)) return args[i+1];
  return def;
}

