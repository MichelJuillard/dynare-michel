#include "command_line_VAR.h"
#include "VARio.h"
#include "switchio.h"
#include "dw_error.h"
#include "dw_parse_cmd.h"
#include "dw_ascii.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/*
  Allocates memory for filename.  Assumes that fmt is of the form

     %s*%s*

  where the first %s will be filled with dir and the second will be
  filled with tag.  If either dir or tag is a null pointer, then the
  the empty string will be used.  The calling routine must free the
  returned pointer.
*/
char* CreateFilenameFromTag(char *fmt, char *tag, char *dir)
{
  char *filename;
  if (!tag) tag="";
  if (!dir) dir="";
  sprintf(filename=(char*)malloc(strlen(dir) + strlen(fmt) + strlen(tag) - 3),fmt,dir,tag);
  return filename;
}

/*
   Create a full path name by appending a "/" if necessary.  The 
   returned pathname must be freed by he calling routine.
*/
char* CreatePath(char *path)
{
#define DIR_DELIMITER '\\'
  char *fullpath;
  int n;
  if (!path) path="";
  n=(int)strlen(path);
  if (path[0] && path[n-1] != DIR_DELIMITER)
    {
      memcpy(fullpath=(char*)malloc(n+2),path,n);
      fullpath[n]=DIR_DELIMITER;
      fullpath[n+1]='\0';
    }
  else
    fullpath=dw_DuplicateString(path);
  return fullpath;
#undef DIR_DELIMITER
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
static char *help_options[]={"-di","-do","-fs","-fp","-ph","-pho","-MLE","-ft","-fto",(char*)NULL};
static char *help_messages[]=
{
  "-di <directory>",
  "If this argument exists, then all input files are in specified directory.",
  "-do <directory>",
  "If this argument exists, then all output files are in specified directory.",
  "-fs <filename>",
  "If this argument exists, then the specification file name is <filename>.",
  "-fp <filename>",
  "If this argument exists, then the parameter file name is <filename>.",
  "-ph <header>",
  "If this argument exists, then the parameter header is <header>.  The default value is \"Posterior mode: \", unless -MLE is in the command line, in which case it is \"MLE: \".",
  "-pho <header>",
  "If this argument exists, then the parameter header used for output is <header>.  The default value is -ph <header>.",
  "-MLE",
  "If this augument exists, them \"MLE: \" is the default value for -ph.",
  "-ft <tag>",
  "The input file tag.  Used to create input filenames if the -fs or -fp options are not present.",
  "-fto <tag>",
  "The output file tag.  Used to create output filenames.  The default value is -ft <tag>.",
  (char*)NULL,
  (char*)NULL
};

static int Match(char *option, char **list, int step)
{
  int n=0, i=0;
  while (option[n] && !isspace(option[n])) n++;
  if (n > 0)
  while (list[i])
    if (memcmp(option,list[i],n))
      i+=step;
    else
      return i;
  return -1;
}

static void PrintMessage(FILE *f, char *msg)
{
#define LL 76
  int k;
  char line[LL];
  while (1)
    {
      while (*msg && isspace(*msg)) msg++;
      if (!(*msg)) return;
      strncpy(line,msg,LL);
      k=LL-1;
      if (!line[k])
        {
          fprintf(f,"    %s\n",line);
          return;
        }
      if (isspace(line[k]))
        {
          line[k]='\0';
          msg+=k+1;
        }
      else
        {
          for ( ; k > 0 && !isspace(line[k-1]); k--);
          if (k == 0) k=LL-1;
          line[k]='\0';
          msg+=k;
        }
      fprintf(f,"    %s\n",line);
    }
#undef LL
}

void PrintHelpMessages(FILE *f, char **include, char **additional)
{

  int i, j;

  if (!f) return;

  if (include)
    {
      i=0;
      while (include[i])
        if ((j=Match(include[i++],help_messages,2)) != -1)
          {
            fprintf(f,"  %s\n",help_messages[j]);
            PrintMessage(f,help_messages[j+1]);
          }
    }
  else
    {
      i=0;
      while (help_messages[i])
        {
          fprintf(f,"  %s\n",help_messages[i++]);
          PrintMessage(f,help_messages[i++]);
        }
    }

  if (additional)
    {
      i=0;
      while (additional[i])
        {
          fprintf(f,"  %s\n",additional[i++]);
          PrintMessage(f,additional[i++]);
       }
    }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

TVARCommandLine* Create_VARCommandLine(void)
{
  TVARCommandLine *cmd=(TVARCommandLine*)malloc(sizeof(TVARCommandLine));
  if (cmd)
    {
       cmd->out_directory=(char*)NULL;
       cmd->in_directory=(char*)NULL;
       cmd->in_tag=(char*)NULL;
       cmd->out_tag=(char*)NULL;
       cmd->out_header=(char*)NULL;

       cmd->specification_filename=(char*)NULL;
       cmd->parameters_filename=(char*)NULL;
       cmd->parameters_header=(char*)NULL;

       cmd->specification_filename_actual=(char*)NULL;
       cmd->parameters_filename_actual=(char*)NULL;
       cmd->parameters_header_actual=(char*)NULL;
    }

  return cmd;
}

void Free_VARCommandLine(TVARCommandLine *cmd)
{
  if (cmd)
    {
      if (cmd->out_directory) free(cmd->out_directory);
      if (cmd->in_directory) free(cmd->in_directory);
      if (cmd->in_tag) free(cmd->in_tag);
      if (cmd->out_tag) free(cmd->out_tag);
      if (cmd->out_header) free(cmd->out_header);
      if (cmd->specification_filename) free(cmd->specification_filename);
      if (cmd->parameters_filename) free(cmd->parameters_filename);
      if (cmd->parameters_header) free(cmd->parameters_header);   
      if (cmd->specification_filename_actual) free(cmd->specification_filename_actual);
      if (cmd->parameters_filename_actual) free(cmd->parameters_filename_actual);
      if (cmd->parameters_header_actual) free(cmd->parameters_header_actual);   
    }
}

/*
  -di <directory>
      If this argument exists, then all input files are in specified directory.

   -do <directory>
      If this argument exists, then all output files are in specified directory.

   -fs <filename>
      If this argument exists, then the specification file name is <filename>.  
      The argument -fs takes precedence over -ft.

   -fp <filename>
      If this argument exists, then the parameter file name is <filename>.  The 
      default value is the filename associated with the argument -fs or -ft.

   -ph <header>
      If this argument exists, then the header for the parameters is <header>.  
      The default value is "MLE: " if -MLE is in the command line and 
      "Posterior mode: " otherwise.

   -pho <header>
      If this argument exists, then the parameter header used for output is 
      <header>.  The default value is -ph <header>.

   -MLE
      If this augument exists, default value for the parameters header is "MLE: ".
      This assumes that the estimate file was produced via a maximum likelihood
      estimate.

   -ft <tag>
      The input file tag.  Used to create input filenames if the -fs or
      -fp options are not present.
    
   -fto <tag>
      The output file tag.  Used to create output filenames.  The default
      value is -ft <tag>.
*/
TVARCommandLine* Base_VARCommandLine(int nargs, char **args, TVARCommandLine *cmd)
{
  if (!cmd && !(cmd=Create_VARCommandLine())) return (TVARCommandLine*)NULL;

  // input directory
  cmd->in_directory=CreatePath(dw_ParseString_String(nargs,args,"di",""));

  // output directory
  cmd->out_directory=CreatePath(dw_ParseString_String(nargs,args,"do",""));
  
  // Specification file
  cmd->specification_filename=dw_DuplicateString(dw_ParseString_String(nargs,args,"fs",(char*)NULL));

  // Parameters file
  cmd->parameters_filename=dw_DuplicateString(dw_ParseString_String(nargs,args,"fp",(char*)NULL));

  // Parameter header
  cmd->parameters_header=dw_DuplicateString((dw_FindArgument_String(nargs,args,"MLE") == -1)
     ? dw_ParseString_String(nargs,args,"ph","Posterior mode: ")
     : dw_ParseString_String(nargs,args,"ph","MLE: "));
						 
  // Output parameters header
  cmd->out_header=dw_DuplicateString(dw_ParseString_String(nargs,args,"ph",cmd->parameters_header));

  // Input file tag
  cmd->in_tag=dw_DuplicateString(dw_ParseString_String(nargs,args,"ft",(char*)NULL));

  // Output file tag
  cmd->out_tag=dw_DuplicateString(dw_ParseString_String(nargs,args,"fto",cmd->in_tag));

  return cmd;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
TStateModel* CreateTStateModelFromEstimateFinal(int nargs, char **args, TVARCommandLine **p_cmd)
{
  TStateModel *model=(TStateModel*)NULL;
  char *filename;
  TVARCommandLine *cmd;

  if (!(*p_cmd) && !(*p_cmd=Base_VARCommandLine(nargs,args,*p_cmd))) return (TStateModel*)NULL;
  cmd=*p_cmd;

  if (cmd->specification_filename)
    {
      filename=CreateFilenameFromTag("%s%s",cmd->specification_filename,cmd->in_directory);
      if (!(model=Read_VAR_Specification((FILE*)NULL,filename)))
        {
          free(filename);
          return (TStateModel*)NULL;
        }
    }
  else
    if (cmd->in_tag)
      {
        filename=CreateFilenameFromTag("%sest_final_%s.dat",cmd->in_tag,cmd->in_directory);
        if (!(model=Read_VAR_Specification((FILE*)NULL,filename)))
          {
            free(filename);
            return (TStateModel*)NULL;
          }
      }
    else
      return (TStateModel*)NULL;

  if (cmd->specification_filename_actual) free(cmd->specification_filename_actual);
  cmd->specification_filename_actual=filename;

  filename=(cmd->parameters_filename) 
    ? CreateFilenameFromTag("%s%s",cmd->parameters_filename,cmd->in_directory)
    : dw_DuplicateString(cmd->specification_filename_actual);

  if (!ReadTransitionMatrices((FILE*)NULL,filename,cmd->parameters_header,model) 
              || !Read_VAR_Parameters((FILE*)NULL,filename,cmd->parameters_header,model))
    {
      free(filename);
      FreeStateModel(model);
      return (TStateModel*)NULL;
    }

  if (cmd->parameters_filename_actual) free(cmd->parameters_filename_actual);
  cmd->parameters_filename_actual=filename;
  if (cmd->parameters_header_actual) free(cmd->parameters_header_actual);
  cmd->parameters_header_actual=dw_DuplicateString(cmd->parameters_header);

  return model;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*
   Attempts to get the parameters from the last iteration in the intermediate 
   file.  Returns one and sets cmd->parameters_file_actual, 
   cmd->parameters_header_actual and loads parameters upon success.  Returns 
   zero upon failure.
*/
int GetLastIteration(TStateModel *model, TVARCommandLine *cmd)
{
  char *filename, *header, *fmt="Iteration %d: ";
  int rtrn=0, cont, terminal_errors, i, j, k=1;
  FILE *f_in;

  filename=CreateFilenameFromTag("%sest_intermediate_%s.dat",cmd->in_tag,cmd->in_directory);
  if (!(f_in=fopen(filename,"rt")))
    {
      free(filename);
      return 0;
    }

  terminal_errors=dw_SetTerminalErrors(dw_GetTerminalErrors() & (~USER_ERR));

  do    
    {
      for (j=10, i=1; k >= j; j*=10, i++);
      sprintf(header=(char*)malloc(strlen(fmt) + i - 1),fmt,k);
      if (ReadTransitionMatrices(f_in,(char*)NULL,header,model) 
                          && Read_VAR_Parameters(f_in,(char*)NULL,header,model))
        {
          cont=1;
          free(header);
          k++;
        }
      else
        cont=0;
    }
 while (cont);

 if (k > 1)
   {
     k--;
     for (j=10, i=1; k >= j; j*=10, i++);
     sprintf(header=(char*)malloc(strlen(fmt) + i - 1),fmt,k);
     if (ReadTransitionMatrices(f_in,(char*)NULL,header,model) && Read_VAR_Parameters(f_in,(char*)NULL,header,model))
       {
         if (cmd->parameters_filename_actual) free(cmd->parameters_filename_actual);
         cmd->parameters_filename_actual=filename;
         if (cmd->parameters_header_actual) free(cmd->parameters_header_actual);
         cmd->parameters_header_actual=header;
         dw_SetTerminalErrors(terminal_errors);
         return 1;
       }
   }
 else
   {
     header="Initial: ";
     if (ReadTransitionMatrices(f_in,(char*)NULL,header,model)  && Read_VAR_Parameters(f_in,(char*)NULL,header,model))
       {
         if (cmd->parameters_filename_actual) free(cmd->parameters_filename_actual);
         cmd->parameters_filename_actual=filename;
         if (cmd->parameters_header_actual) free(cmd->parameters_header_actual);
         cmd->parameters_header_actual=dw_DuplicateString(header);
         dw_SetTerminalErrors(terminal_errors);
         return 1;
       }
   }

 free(filename);
 dw_SetTerminalErrors(terminal_errors);
 return 0;
}

/*
   Attempt to set up model from command line. 

   -ft <filename tag>
      If this argument exists, then the following is attempted:

         specification file name:  est_final_<tag>.dat
         init/restart file name:   est_final_<tag>.dat with header="Posterior mode: "

         specification file name:  init_<tag>.dat
         init/restart file name:   est_intermediate_<tag>.dat with header="Iteration %d: "

         (not yet implemented)
         specification file name:  init_<tag>.dat
         init/restart file name:   est_csminwel_<tag>.dat  

         specification file name:  init_<tag>.dat
         init/restart file name:   init_<tag>.dat with header="Initial: "
    
   Returns valid pointer to a TStateModel upon success and null upon failure.
*/
TStateModel* CreateTStateModelForEstimate(int nargs, char **args, TVARCommandLine **p_cmd)
{
  TStateModel *model;
  char *filename, *header;
  int terminal_errors;
  TVARCommandLine *cmd;

  terminal_errors=dw_SetTerminalErrors(dw_GetTerminalErrors() & (~USER_ERR));

  if (!(*p_cmd) && !(*p_cmd=Base_VARCommandLine(nargs,args,*p_cmd))) return (TStateModel*)NULL;
  cmd=*p_cmd;

  if (cmd->specification_filename)
    {
      filename=CreateFilenameFromTag("%s%s",cmd->specification_filename,cmd->in_directory);
      if (!(model=Read_VAR_Specification((FILE*)NULL,filename)))
        {
          free(filename);
          dw_SetTerminalErrors(terminal_errors);
          return (TStateModel*)NULL;
        }
    }
  else
    if (cmd->in_tag)
      {
        filename=CreateFilenameFromTag("%sest_final_%s.dat",cmd->in_tag,cmd->in_directory);
        if (!(model=Read_VAR_Specification((FILE*)NULL,filename)))
          {
            free(filename);
            filename=CreateFilenameFromTag("%sinit_%s.dat",cmd->in_tag,cmd->in_directory);
            if (!(model=Read_VAR_Specification((FILE*)NULL,filename)))
              {
                free(filename);
                dw_SetTerminalErrors(terminal_errors);
                return (TStateModel*)NULL;
              }
          }
      }
    else
      {
        dw_SetTerminalErrors(terminal_errors);
        return (TStateModel*)NULL;
      }

  if (cmd->specification_filename_actual) free(cmd->specification_filename_actual);
  cmd->specification_filename_actual=filename;
 
  if (cmd->parameters_filename)
    {
      header=cmd->parameters_header;
      filename=CreateFilenameFromTag("%s%s",cmd->parameters_filename,cmd->in_directory);
      if (!ReadTransitionMatrices((FILE*)NULL,filename,header,model) 
		  || !Read_VAR_Parameters((FILE*)NULL,filename,header,model))
        {
          free(filename);
          FreeStateModel(model);
          dw_SetTerminalErrors(terminal_errors);
          return (TStateModel*)NULL;
        }
    }
  else
    if (cmd->specification_filename)
      {
        header=cmd->parameters_header;
        filename=CreateFilenameFromTag("%s%s",cmd->specification_filename,cmd->in_directory);
        if (!ReadTransitionMatrices((FILE*)NULL,filename,header,model) 
		    || !Read_VAR_Parameters((FILE*)NULL,filename,header,model))
          {
            free(filename);
            FreeStateModel(model);
            dw_SetTerminalErrors(terminal_errors);
            return (TStateModel*)NULL;
          }
      }
    else
      if (cmd->in_tag)
        {
          header=cmd->parameters_header;
          filename=filename=CreateFilenameFromTag("%sest_final_%s.dat",cmd->in_tag,cmd->in_directory);
          if (!ReadTransitionMatrices((FILE*)NULL,filename,header,model)
                      || !Read_VAR_Parameters((FILE*)NULL,filename,header,model))
            {
              free(filename);
              if (GetLastIteration(model,cmd))
                {
                  dw_SetTerminalErrors(terminal_errors);
                  return model;
                }
              else
                {                                   
                  header="Initial: ";
                  filename=filename=CreateFilenameFromTag("%sinit_%s.dat",cmd->in_tag,cmd->in_directory);
                  if (!ReadTransitionMatrices((FILE*)NULL,filename,header,model)
                              || !Read_VAR_Parameters((FILE*)NULL,filename,header,model))
                    {
                      FreeStateModel(model);
                      dw_SetTerminalErrors(terminal_errors);
                      return (TStateModel*)NULL;
                    }
                }
            }
        }

  if (cmd->parameters_filename_actual) free(cmd->parameters_filename_actual);
  cmd->parameters_filename_actual=filename;
  if (cmd->parameters_header_actual) free(cmd->parameters_header_actual);
  cmd->parameters_header_actual=dw_DuplicateString(header);

  dw_SetTerminalErrors(terminal_errors);
  return model;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
