
#include "dw_ascii.h"
#include "dw_array.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>

#include "modify_for_mex.h"

/*
   Attempts to open filename for reading.  Returns pointer to file upon success 
   and prints error message and exits upon failure.  The file must exist.
*/
FILE *dw_OpenTextFile(char *filename)
{
  FILE *f=fopen(filename,"rt");
  if (!f)
    {
      printf("Unable to open %s\n",filename);
      exit(0);
    }
  return (f);
}

/*
   Attempts to create filename for writing.  Returns pointer to file upon success 
   and prints error message and exits upon failure.  If the file exists, it is
   overwritten.
*/
FILE *dw_CreateTextFile(char *filename)
{
  FILE *f=fopen(filename,"wt");
  if (!f)
    {
      printf("Unable to open %s\n",filename);
      exit(0);
    }
  return (f);
}

/*
   Attempts to create filename for writing.  Returns pointer to file upon success 
   and prints error message and exits upon failure.  The file is created if it 
   does not exist and is opened with the file pointer positioned at the end of
   file if it does exist.
*/
FILE *dw_AppendTextFile(char *filename)
{
  FILE *f=fopen(filename,"at");
  if (!f)
    {
      printf("Unable to open %s\n",filename);
      exit(0);
    }
  return (f);
}

/*
   Assumes:
     f      : valid file pointer
     buffer : pointer to character or null pointer
     n      : pointer to integer containing the length of buffer

   Returns:
     Pointer to null terminated string containing the characters from the file up
     to and including the terminating new line character.  A null pointer return
     indicates that there was a memory error or no characters to be read.  Call
     dw_GetError() to determine if a error occured.   

   Results:
     Reads line, beginning at current position from file f  Returns a pointer to 
     the buffer containing the file and resets *n if necessary.  The if the 
     passed buffer is null or is not large enough to contain the line, buffer is 
     freed and a new buffer is allocated.  Because of this, the passed buffer 
     must either null or allocated with malloc(), realloc(), or calloc() and the 
     calling routine is responsible for eventually freeing the memory if the 
     return value is not null.

   Notes:
     If buffer is null, then value pointed to by the pointer n is not used.
*/
#define SIZE_INCREMENT 1024
char* dw_ReadLine(FILE *f, char *buffer, int *n)
{
  char *ptr, *nbuffer;
  int i, k=0;
  if (!buffer && !(buffer=(char*)malloc(*n=SIZE_INCREMENT)))
    {
      *n=0;
      return (char*)NULL;
    }
  ptr=buffer;
  while (fgets(ptr,*n-k,f))
    if (ptr[(i=(int)strlen(ptr))-1] == '\n')
      return buffer;
    else
      if (!(nbuffer=(char*)realloc(buffer,*n+=SIZE_INCREMENT)))
	{
	  free(buffer);
	  *n=0;
	  return (char*)NULL;
	}
      else
	ptr=(buffer=nbuffer) + (k+=i);
  if (ptr != buffer)
    return buffer;
  else
    {
      free(buffer);
      *n=0;
      return (char*)NULL;
    }
}
#undef SIZE_INCREMENT


char** dw_ParseDelimitedString(char *buffer, char delimiter, int flag)
{
  struct StringList
  {
    struct StringList *next;
    char *string;
    int length;
  } *head, *ptr;
  int k=0, n, m;
  char **v;
  if (!buffer) return (char**)NULL;
  for (head=ptr=(struct StringList*)NULL; *buffer; buffer+=buffer[n] ? n+1 : n)
    {
      if (flag & STRIP_LEADING_WHITESPACE)
	while (*buffer && (*buffer != delimiter) && isspace(*buffer)) buffer++;
      for (n=0; buffer[n] && (buffer[n] != delimiter); n++);
      if (flag & STRIP_TRAILING_WHITESPACE)
        for (m=n-1; (m >= 0) && isspace(buffer[m]); m--);
      else
        m=n-1;
      if ((m >= 0) || !(flag & REMOVE_EMPTY_FIELDS))
        {
          ptr=(struct StringList*)malloc(sizeof(struct StringList));
          ptr->string=buffer;
          ptr->length=m+1;
          ptr->next=head;
          head=ptr;
          k++;
        }
    }
  v=dw_CreateArray_string(k);
  while (--k >= 0)
    {
      v[k]=(char*)malloc(head->length+1);
      if (head->length > 0) memcpy(v[k],head->string,head->length);
      v[k][head->length]='\0';
      ptr=head;
      head=head->next;
      free(ptr);
    }
  return v;
}

/*
   Assumes
     f:  valid file pointer
     delimiter:  field deliniter.
     flag:  one of the values defined in dw_ascii.h

   Returns
     One-dimensional string array of the delimited fields of the current line of 
     the file f or a null pointer.

   Notes
     The file is read starting from the current file position.  If the file 
     contains no fields or there is a memory error, then a null pointer is 
     returned.  The delimiter character defines the fields in each row and the 
     new line character defines the rows. 
*/
char** dw_ReadDelimitedLine(FILE *f, char delimiter, int flag)
{
  int n=0;
  char **v=(char**)NULL, *buffer=dw_ReadLine(f,(char*)NULL,&n);
  if (buffer)
    {
      v=dw_ParseDelimitedString(buffer,delimiter,flag);
      free(buffer);
    }
  return v;
}

/*
   Assumes
     f:  valid file pointer or null pointer.
     filename:  pointer to null terminated string or null pointer.
     delimiter:  field deliniter.
     flag:  one of the values defined in dw_ascii.h

   Returns
     Two-dimensional string array of the deliminted fields of f or a null 
     pointer.

   Notes
     One of f and filename should be non-null.  If f is non-null, the file is 
     read starting from the current file position.  If f is null, an attempt is
     made to open the file.  If successful, the file is read from the beginning.
     If the file does not exist or contains no fields, then a null pointer is 
     returned.  The delimiter character defines the fields in each row and the 
     new line character defines the rows. 

*/
char*** dw_ReadDelimitedFile(FILE *f, char* filename, char delimiter, int flag)
{
  struct LineList
    {
      struct LineList *next;
      char **line;
    } *head=(struct LineList*)NULL, *ptr;
  int n=0;
  char **v, ***M=(char***)NULL, *buffer=(char*)NULL;
  FILE *f_in=f ? f : fopen(filename,"rt");
  if (f_in)
    {
      while (buffer=dw_ReadLine(f_in,buffer,&n))
	if (v=dw_ParseDelimitedString(buffer,delimiter,flag))
          {
	    ptr=(struct LineList*)malloc(sizeof(struct LineList));
	    ptr->line=v;
            ptr->next=head;
            head=ptr;
            n++;
          }
      if (!f) fclose(f_in);
      if (n > 0)
        {
          M=(char***)dw_CreateArray_array(n);
          while (--n >= 0)
            {
              M[n]=head->line;
              ptr=head;
              head=head->next;
              free(ptr);
            }
        }
    }
  return M;
}

int dw_PrintDelimitedArray(FILE *f, void* array, char delimiter)
{
  char format[4];
  format[0]='%';
  format[1]='s';
  format[2]=delimiter;
  format[3]='\0';
  return dw_PrintArray(f,array,format);
}


/*
   Assumes:
     f         : valid file pointer
     delimiter : field terminator
     terminal  : line terminator
     flag      : determine how characters are processed
     buffer    : pointer to pointer to character or null pointer
     n         : pointer to integer containing the length of buffer

   Returns:
     0 : memory error occured
     1 : field read, terminated by delimiter
     2 : field read, terminated by terminal
     3 : field read, terminated by EOF

   Results:
     If necessary, memory ia reallocated.  The length of this reallocated memory 
     is stored in n.  It is the calling routines responsibility to free the
     memory pointed to by *buffer.  

   Notes:
     flag values
       ALLOW_QUOTED_TEXT
         If set the delimiter and terminal characters do not stop processing when
         encountered between quotes.  To produce a quote in quoted text, use two
         consectutive quotes.  Outside quoted text, a quote always begins quoted
         text.

       PRINTABLE_ONLY_IN_QUOTES

       PRINTABLE_ONLY

       STRIP_LEADING_WHITESPACE

       STRIP_TRAILING_WHITESPACE

       STRIP_WHITESPACE


*/
//#define INCREMENT 1024
//int dw_ReadDelimitedField(FILE *f, int delimiter, int terminal, int flag, char **buffer, int *n)
//{
/*   int ch;        // next character read */
/*   int k=0;       // position to store next char, always less than *n */
/*   int quoted=0; */
/*   int leading=(flag & STRIP_LEADING_WHITESPACE) ? 1 : 0; */
/*   char *ptr; */

/*   ch=fgetc(f); */

/*   while (ch != EOF) */
/*     { */
/*       //=== reallocate memory if necessary */
/*       if (k+1 > *n) */
/* 	if (!(ptr=(char*)realloc(buffer,*n+=INCREMENT))) */
/* 	  { */
/* 	    *n-=INCREMENT; */
/* 	    return 0; */
/* 	  } */
/* 	else */
/* 	  buffer=ptr; */

/*       //=== process character */
/*       if (quoted) */
/* 	{ */
/* 	  if (ch == '"') */
/* 	    if ((ch=fgets(f)) != '"')  */
/* 	      { */
/* 		quoted=0; */
/* 		continue; */
/* 	      } */
/* 	  if (!(flag & PRINTABLE_ONLY_IN_QUOTES) || isprint(ch))  */
/* 	    buffer[k++]=ch; */
/* 	} */
/*       else */
/* 	if ((ch == delimiter) || (ch == terminal))  */
/* 	  break; */
/* 	else */
/* 	  if ((ch == '"') && (flag & ALLOW_QUOTED_TEXT)) */
/* 	    quoted=1; */
/* 	  else */
/* 	    if (!(flag & PRINTABLE_ONLY) || isprint(ch)) */
/* 	      { */
/* 		if ((ch == "\r") && (terminal == '\n')) */
/* 		  { */
/* 		    if ((ch=fgetc(f)) == '\n') break; */
/* 		    if (!leading) buffer[k++]='\r'; */
/* 		    continue; */
/* 		  } */
/* 		if (leading) */
/* 		  if (isspace(ch))  */
/* 		    { */
/* 		      ch=fgetc(f); */
/* 		      continue; */
/* 		    } */
/* 		  else */
/* 		    leading=0; */
/* 		buffer[k++]=ch; */
/* 	      } */

/*       ch=fgets(f); */
/*     } */

/*   buffer[k]='\0'; */

/*   return (ch == EOF) ? 3 : (ch == terminal) ? 2 : 1;  */
//}
//#undef INCREMENT

/*
   Returns 1 if the null terminated string id is found at the beginning of a line
   in the file and 0 otherwise.  The file pointer is set to the line immediately 
   after the line containing id.  The search starts at the current position of
   the file.  If id is not found, then the file is rewound and the search is 
   continued until the initial file position is passed.
*/
int dw_SetFilePosition(FILE *f, char *id)
{
  char *buffer=(char*)NULL;
  int m, n, pos;
  if ((n=(int)strlen(id)) > 0)
    {
      pos=ftell(f);
      while (buffer=dw_ReadLine(f,buffer,&m))
	if (!memcmp(buffer,id,n)) 
	  {
	    free(buffer);
	    return 1;
	  }
      if (pos > 0)
	{
	  rewind(f);
	  while ((ftell(f) < pos) && (buffer=dw_ReadLine(f,buffer,&m)))
	    if (!memcmp(buffer,id,n))
	      {
		free(buffer);
		return 1;
	      }
	  if (buffer) free(buffer);
	}
    }
  return 0;
}

/*
   Returns 1 if the null terminated string id is found at the beginning of a line
   in the file and 0 otherwise.  The file pointer is set to the line immediately 
   after the line containing id.  Compares a maximum of 1023 characters of id.  
   The file is not rewound so that the search starts at the current position.
*/
int dw_SetFilePositionNoRewind(FILE *f, char *id)
{
 char buffer[1024], ch;
 int n=(int)strlen(id);
 if (n > 1023) n=1023;
 while (fgets(buffer,1024,f))
  {
   if (buffer[strlen(buffer)-1] != '\n')
    do 
     ch=fgetc(f);
    while ((ch != '\n') && (ch != EOF));
   if (!memcmp(buffer,id,n)) return 1;
  }
 return 0;
}


int dw_SetFilePositionBySection(FILE *f, int n, ...)
{
  char *arg;
  int i;
  va_list ap;
  rewind(f);
  va_start(ap,n);
  for (i=0; i  < n; i++)
    if (!(arg=va_arg(ap,char*)) || !dw_SetFilePositionNoRewind(f,arg))
      {
	va_end(ap);
	return 0;
      }
  va_end(ap);
  return 1;
}

char* dw_DuplicateString(char *buffer)
{
  char *rtrn=(char*)NULL;
  if (buffer && (rtrn=(char*)malloc(strlen(buffer)+1))) strcpy(rtrn,buffer);
  return rtrn;
}
