
#include "dw_array.h"
#include "dw_error.h"

#include <stdlib.h>
#include <string.h>

#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif
#include <stdarg.h>

#include "modify_for_mex.h"

/*  //================================== Macros ===================================//   ansi-c*/
#define dw_ElementSizeA(a)                (dw_SpecsA(a)->size)
#define dw_GetOffsetA(a)                  (dw_SpecsA(a)->offset)
#define dw_IsSameTypeA(a1,a2)             (!memcmp(dw_SpecsA(a1),dw_SpecsA(a2),sizeof(TElementSpecification)))
#define dw_IsPointerA(a)                  (dw_SpecsA(a)->flag & dw_ARRAY_POINTER)
#define dw_UseMemcpyA(a)                  (dw_SpecsA(a)->flag & dw_ARRAY_USE_MEMCPY)
#define dw_DeleteSpecsA(a)                (dw_SpecsA(a)->flag & dw_ARRAY_DELETE_SPECS)
#define dw_GetDestructorA(a)              (dw_SpecsA(a)->destructor)
#define dw_GetDefaultConstructorA(s)      (dw_SpecsA(a)->default_constructor)
#define dw_GetPointerCopyConstructorA(a)  (dw_SpecsA(a)->pointer_copy_constructor)
#define dw_GetStaticCopyConstructorA(a)   (dw_SpecsA(a)->static_copy_constructor)
#define dw_GetPrintRoutineA(a)            (dw_SpecsA(a)->print_routine)
#define dw_GetReadRoutineA(a)             (dw_SpecsA(a)->read_routine)


/*******************************************************************************/
/********************** C-style multi-dimensional arrays ***********************/
/*******************************************************************************/
/*
   Frees a C-style multi-dimensional array.  The pointer a must point to a valid
   array created via a call to dw_CreateArray() or be a null pointer.
*/
void dw_FreeArray(void* a)
{
  int i, size, offset;
  void (*Destructor)(void*);
  if (a)
    {
      if (Destructor=dw_GetDestructorA(a))
    if (dw_IsPointerA(a))
      for (i=dw_DimA(a)-1; i >= 0; i--)
        Destructor(((void**)a)[i]);
    else
      for (i=(size=dw_ElementSizeA(a))*(dw_DimA(a)-1); i >= 0; i-=size)
        Destructor((void*)(((char*)a) + i));
      offset=dw_GetOffsetA(a);
      if (dw_DeleteSpecsA(a)) free(dw_SpecsA(a));
      free((void*)(((char*)a) - offset));
    }
}

/*
   Assumes:
     specs:  Pointer to a valid TElementSpecification structure.
     dim:  Positive integer

   Returns:
     A pointer to a valid array of lenth dim upon success and a null pointer upon
     failure.

   Notes:
     The return value should be type cast to the appropriate pointer type.
*/
void* dw_CreateArray(TElementSpecification *specs, int dim)
{
  void *a=(void*)NULL;
  int i;
  if (dim <= 0)
    dw_Error(ARG_ERR);
  else
    if (!(a=swzMalloc(dim*specs->size + specs->offset)))
      dw_Error(MEM_ERR);
    else
      {
    a=(void*)(((char*)a)+specs->offset);
    dw_DimA(a)=dim;
    dw_SpecsA(a)=specs;
    if (specs->default_constructor)
      for (i=(specs->size)*(dim-1); i >= 0; i-=specs->size)
        specs->default_constructor((void*)(((char*)a) + i));
      }
  return a;
}

/*
   Assumes:
     specs:  Pointer to a valid TElementSpecification structure.
     depth:  Positive integer
     dim:  Array of positive integers of length at least depth

   Returns:
     A pointer to a valid multidimensiona array.  The dimensions of the array are
     determined by depth and dim.

   Notes:
     The return value should be type cast to the appropriate pointer type.
*/
void* dw_CreateMultidimensionalArray(TElementSpecification *specs, int depth, int *dim)
{
  int i;
  void *a;
  if (depth == 1) return dw_CreateArray(specs,dim[0]);
  if (a=dw_CreateArray_array(dim[0]))
    for (i=dim[0]-1; i >= 0; i--)
      if (!(((void**)a)[i]=dw_CreateMultidimensionalArray(specs,depth-1,dim+1)))
    {
      dw_FreeArray(a);
      return (void*)NULL;
    }
  return a;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/**************************** Default Constructors *****************************/
/*******************************************************************************/
void DefaultPointerConstructor(void *element)
{
  *((void**)element)=(void*)NULL;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/******************************* Print Functions *******************************/
/*******************************************************************************/
int dw_PrintArray(FILE *f, void *a, char *format)
{
  int i, size;
  int (*PrintRoutine)(FILE*, void*, char*);
  if (f && a)
    if (PrintRoutine=dw_GetPrintRoutineA(a))
      {
    if (dw_IsPointerA(a))
      for (i=0; i < dw_DimA(a); i++)
        { if (!PrintRoutine(f,((void**)a)[i],format)) return 0; }
    else
      for (size=dw_ElementSizeA(a), i=0; i < dw_DimA(a); i++)
        { if (!PrintRoutine(f,(void*)(((char*)a) + i*size),format)) return 0; }

        if(f==stdout)
          printf("\n");
        else
          fprintf(f,"\n");
    return 1;
      }
  return 0;
}

static int dw_PrintInt(FILE* f, void* element, char *format)
{
  if(f==stdout)
    {
      swz_fprintf_stdout(format ? format : "%d ",*((int*)element));
      return 1;
    }
  else
    return (fprintf(f,format ? format : "%d ",*((int*)element)) < 0) ? 0 : 1;
}

static int dw_PrintDouble(FILE* f, void* element, char *format)
{
  if(f==stdout)
    {
      swz_fprintf_stdout(format ? format : "%lf ",*((double*)element));
      return 1;
    }
  else
    return (fprintf(f,format ? format : "%lf ",*((double*)element)) < 0) ? 0 : 1;
}

static int dw_PrintFloat(FILE* f, void* element, char *format)
{
  if(f==stdout)
    {
      swz_fprintf_stdout(format ? format : "%f ",*((float*)element));
      return 1;
    }
  else
    return (fprintf(f,format ? format : "%f ",*((float*)element)) < 0) ? 0 : 1;
}

static int dw_PrintChar(FILE* f, void* element, char *format)
{
  if(f==stdout)
    {
      swz_fprintf_stdout(format ? format : "%c ",*((char*)element));
      return 1;
    }
  else
    return (fprintf(f,format ? format : "%c ",*((char*)element)) < 0) ? 0 : 1;
}

static int dw_PrintString(FILE* f, void* element, char *format)
{
  if(f==stdout)
    {
      swz_fprintf_stdout(format ? format : "%s\t",(char*)element);
      return 1;
    }
  else
    return (fprintf(f,format ? format : "%s\t",(char*)element) < 0) ? 0 : 1;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/******************************* Read Functions ********************************/
/*******************************************************************************/
int dw_ReadArray(FILE *f, void *a)
{
  int i, size;
  int (*ReadRoutine)(FILE*, void*);
  if (f && a)
    if (ReadRoutine=dw_GetReadRoutineA(a))
      {
    if (dw_IsPointerA(a))
      for (i=0; i < dw_DimA(a); i++)
        { if (!ReadRoutine(f,((void**)a)[i])) return 0; }
    else
      for (size=dw_ElementSizeA(a), i=0; i < dw_DimA(a); i++)
        { if (!ReadRoutine(f,(void*)(((char*)a) + i*size))) return 0; }
    return 1;
      }
  return 0;
}

static int dw_ReadInt(FILE* f, void* element)
{
  return (fscanf(f," %d ",(int*)element) != 1) ? 0 : 1;
}

static int dw_ReadDouble(FILE* f, void* element)
{
  return (fscanf(f," %lf ",(double*)element) != 1) ? 0 : 1;
}

static int dw_ReadFloat(FILE* f, void* element)
{
return (fscanf(f," %f ",(float*)element) != 1) ? 0 : 1;
}

static int dw_ReadChar(FILE* f, void* element)
{
return (fscanf(f," %c ",(char*)element) != 1) ? 0 : 1;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/****************************** Copy Constructors ******************************/
/*******************************************************************************/
/*
   Assumes
*/
static int FullCopyAttempt(void **d, void *s, void* (*copy)(void*, void*), void (*destructor)(void*))
{
  if (s)
    if (*d)
      {
    if (!copy(*d,s))
      {
        if (destructor) destructor(*d);
        if (!(*d=copy((void*)NULL,s))) return 0;
      }
      }
    else
      {
    if (!(*d=copy((void*)NULL,s))) return 0;
      }
  else
    if (*d)
      {
    if (destructor) destructor(*d);
    *d=(void*)NULL;
      }
  return 1;
}

/*
   Assumes:
     d:  A valid array or null pointer
     s:  A valid array

   Returns:
     Upon success returns a copy of the array s.  If d is null, then the array is
     created.  Upon failure, a null pointer is returned.

   Notes:
     If d is
*/
void* dw_CopyArray(void* d, void* s)
{
  int i, size;
  void* original_d=d;

  if (!s) return (void*)NULL;

  if (s == d) return d;

  if (!d)
    { if (!(d=dw_CreateArray(dw_SpecsA(s),dw_DimA(s)))) return (void*)NULL; }
  else
    { if  ((dw_DimA(s) != dw_DimA(d)) || !dw_IsSameTypeA(d,s)) return (void*)NULL; }

  if (dw_UseMemcpyA(s))
    {
      memcpy(d,s,dw_DimA(s)*dw_ElementSizeA(s));
    }
  else if (dw_GetPointerCopyConstructorA(s))
    {
      for (i=dw_DimA(s)-1; i >= 0; i--)
    if (!FullCopyAttempt(((void**)d)+i,((void**)s)[i],dw_GetPointerCopyConstructorA(s),dw_GetDestructorA(d)))
      {
        if (!original_d) dw_FreeArray(d);
        return (void*)NULL;
      }
    }
  else if (dw_GetStaticCopyConstructorA(s))
    {
      for (i=(size=dw_ElementSizeA(s))*(dw_DimA(s)-1); i >= 0; i-=size)
    if (!dw_GetStaticCopyConstructorA(s)((void*)(((char*)d) + i),(void*)(((char*)s) + i)))
      {
        if (!original_d) dw_FreeArray(d);
        return (void*)NULL;
      }
    }
  else
    {
      if (!original_d) dw_FreeArray(d);
      return (void*)NULL;
    }

  return d;
}

/*
    Assumes
      Both d and s are valid pointers and both *d and *s are either null or a
      null terminated string.  If *d is a null terminated string, then it must
      have been created via a call to swzMalloc(), swzCalloc() or swzRealloc().

    Returns
      Returns one upon success and zero upon failure.

    Results
      If is *s is null, then *d is freed if it is non-null and is then set to
      null.  If *s is null terminated string, then *d is reallocated if more
      memory is required and then *s is copied into *d.

    Notes
      It is critical that this function be called only if the destination string
      was dynamically created via a call to swzMalloc(), swzCalloc() or swzRealloc().  If
      this is not the case, then servere memory problems can result.
*/
static int dw_CopyString(void *d, void *s)
{
  char* dest;
  if (*((char**)s))
    if (dest=swzRealloc(*((char**)d),strlen(*((char**)s))+1))
      strcpy(*((char**)d)=dest,*((char**)s));
    else
      return 0;
  else
    if (*((char**)d))
      {
    free(*((char**)d));
    *((char**)d)=(char*)NULL;
      }
  return 1;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/********* Multidimensional Arrays Create Via Variable Argument Lists **********/
/*******************************************************************************/
/*
   Assumes:
     specs:  Pointer to a valid TElementSpecification structure.
     depth:  Positive integer

   Returns:
     A pointer to a valid multidimensiona array.  The dimensions of the array are
     determined by depth the variable list of arguments.

   Notes:
     The return value should be type cast to the appropriate pointer type.  The
     variable list of arguments must be at least of length depth and consist of
     positive integers.
*/
void* dw_CreateMultidimensionalArrayList(TElementSpecification *specs, int depth, ...)
{
  va_list ap;
  int i, *dim;
  void *a=(void*)NULL;
  if (dim=(int*)swzMalloc(depth*sizeof(int)))
    {
      va_start(ap,depth);
      for (i=0; i < depth; i++) dim[i]=va_arg(ap,int);
      va_end(ap);
      a=dw_CreateMultidimensionalArray(specs,depth,dim);
      free(dim);
    }
  return a;
}

void* dw_CreateMultidimensionalArrayList_string(int depth, ...)
{
  va_list ap;
  int i, *dim;
  void *a=(void*)NULL;
  if (dim=(int*)swzMalloc(depth*sizeof(int)))
    {
      va_start(ap,depth);
      for (i=0; i < depth; i++) dim[i]=va_arg(ap,int);
      va_end(ap);
      a=dw_CreateMultidimensionalArray_string(depth,dim);
      free(dim);
    }
  return a;
}

void* dw_CreateMultidimensionalArrayList_int(int depth, ...)
{
  va_list ap;
  int i, *dim;
  void *a=(void*)NULL;
  if (dim=(int*)swzMalloc(depth*sizeof(int)))
    {
      va_start(ap,depth);
      for (i=0; i < depth; i++) dim[i]=va_arg(ap,int);
      va_end(ap);
      a=dw_CreateMultidimensionalArray_int(depth,dim);
      free(dim);
    }
  return a;
}

void* dw_CreateMultidimensionalArrayList_double(int depth, ...)
{
  va_list ap;
  int i, *dim;
  void *a=(void*)NULL;
  if (dim=(int*)swzMalloc(depth*sizeof(int)))
    {
      va_start(ap,depth);
      for (i=0; i < depth; i++) dim[i]=va_arg(ap,int);
      va_end(ap);
      a=dw_CreateMultidimensionalArray_double(depth,dim);
      free(dim);
    }
  return a;
}

void* dw_CreateMultidimensionalArrayList_float(int depth, ...)
{
  va_list ap;
  int i, *dim;
  void *a=(void*)NULL;
  if (dim=(int*)swzMalloc(depth*sizeof(int)))
    {
      va_start(ap,depth);
      for (i=0; i < depth; i++) dim[i]=va_arg(ap,int);
      va_end(ap);
      a=dw_CreateMultidimensionalArray_float(depth,dim);
      free(dim);
    }
  return a;
}

void* dw_CreateMultidimensionalArrayList_char(int depth, ...)
{
  va_list ap;
  int i, *dim;
  void *a=(void*)NULL;
  if (dim=(int*)swzMalloc(depth*sizeof(int)))
    {
      va_start(ap,depth);
      for (i=0; i < depth; i++) dim[i]=va_arg(ap,int);
      va_end(ap);
      a=dw_CreateMultidimensionalArray_char(depth,dim);
      free(dim);
    }
  return a;
}
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

/*******************************************************************************/
/****************************** Initialize Arrays ******************************/
/*******************************************************************************/
int dw_InitializeArray(void *a, void *x)
{
  int i, size;
  if (a)
    {
      if (dw_IsArrayA(a))
    {
      for (i=dw_DimA(a)-1; i >= 0; i--)
        if (!dw_InitializeArray(((void**)a)[i],x)) return 0;
    }
      else if (dw_UseMemcpyA(a))
    {
      for (size=dw_ElementSizeA(a), i=size*(dw_DimA(a)-1); i >= 0; i-=size)
        memcpy((void*)(((char*)a) + i),x,size);
    }
      else if (dw_GetPointerCopyConstructorA(a))
    {
      for (i=dw_DimA(a)-1; i >= 0; i--)
        if (!FullCopyAttempt(((void**)a)+i,x,dw_GetPointerCopyConstructorA(a),dw_GetDestructorA(a))) return 0;
    }
      else if (dw_GetStaticCopyConstructorA(a))
    {
      for (i=(size=dw_ElementSizeA(a))*(dw_DimA(a)-1); i >= 0; i-=size)
        if (!dw_GetStaticCopyConstructorA(a)((void*)(((char*)a)+i),x)) return 0;
    }
      else
    return 0;
      return 1;
    }
  return 0;
}

int dw_InitializeArray_int(void *a, int x) { return dw_InitializeArray(a,&x); }

int dw_InitializeArray_double(void *a, double x) { return dw_InitializeArray(a,&x); }

int dw_InitializeArray_float(void *a, float x) { return dw_InitializeArray(a,&x); }

int dw_InitializeArray_char(void *a, char x) { return dw_InitializeArray(a,&x); }
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


/*******************************************************************************/
/**************************** TElementSpecification ****************************/
/*******************************************************************************/
TElementSpecification* CreateArraySpecification_pointer(void (*destructor)(void *))
{
  TElementSpecification *specs;
  if (specs=(TElementSpecification*)swzMalloc(sizeof(TElementSpecification)))
    {
      specs->flag=dw_ARRAY_POINTER | dw_ARRAY_DELETE_SPECS;
      specs->size=sizeof(void*);
      specs->offset=sizeof(void*)*((sizeof(int)+sizeof(TElementSpecification*)+sizeof(void*)-1)/sizeof(void*)),
      specs->destructor=destructor;
      specs->default_constructor=NULL;
      specs->pointer_copy_constructor=NULL;
      specs->static_copy_constructor=NULL;
      specs->print_routine=NULL;
      specs->read_routine=NULL;
    }
  return specs;
}

TElementSpecification dw_IntSpecs =
  {
    dw_ARRAY_USE_MEMCPY,
    sizeof(int),
    sizeof(int)*((sizeof(int)+sizeof(TElementSpecification*)+sizeof(int)-1)/sizeof(int)),
    NULL,
    NULL,
    NULL,
    NULL,
    dw_PrintInt,
    dw_ReadInt
  };

TElementSpecification dw_DoubleSpecs =
  {
    dw_ARRAY_USE_MEMCPY,
    sizeof(double),
    sizeof(double)*((sizeof(int)+sizeof(TElementSpecification*)+sizeof(double)-1)/sizeof(double)),
    NULL,
    NULL,
    NULL,
    NULL,
    dw_PrintDouble,
    dw_ReadDouble
  };

TElementSpecification dw_FloatSpecs =
  {
    dw_ARRAY_USE_MEMCPY,
    sizeof(float),
    sizeof(float)*((sizeof(int)+sizeof(TElementSpecification*)+sizeof(float)-1)/sizeof(float)),
    NULL,
    NULL,
    NULL,
    NULL,
    dw_PrintFloat,
    dw_ReadFloat
  };

TElementSpecification dw_CharSpecs =
  {
    dw_ARRAY_USE_MEMCPY,
    sizeof(char),
    sizeof(char)*((sizeof(int)+sizeof(TElementSpecification*)+sizeof(char)-1)/sizeof(char)),
    NULL,
    NULL,
    NULL,
    NULL,
    dw_PrintChar,
    dw_ReadChar
  };

TElementSpecification dw_StringSpecs =
  {
    dw_ARRAY_POINTER,
    sizeof(char*),
    sizeof(char*)*((sizeof(int)+sizeof(TElementSpecification*)+sizeof(char*)-1)/sizeof(char*)),
    free,
    DefaultPointerConstructor,
    NULL,
    dw_CopyString,
    dw_PrintString,
    NULL
  };

TElementSpecification dw_ArraySpecs =
  {
    dw_ARRAY_POINTER | dw_ARRAY_ARRAY,
    sizeof(void*),
    sizeof(void*)*((sizeof(int)+sizeof(TElementSpecification*)+sizeof(void*)-1)/sizeof(void*)),
    dw_FreeArray,
    DefaultPointerConstructor,
    dw_CopyArray,
    NULL,
    dw_PrintArray,
    dw_ReadArray
  };
/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/
