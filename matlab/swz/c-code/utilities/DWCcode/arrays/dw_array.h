
#ifndef __TARRAY__
#define __TARRAY__

/******************************** C-style arrays ********************************

Attempts to implement a C++ class template for multidimensional arrays.  The goal
is to allow access through the bracket operator (a[i_1][i_2]...[i_n]) but provide
mechanisms for creating, destroying, and determining the dimensions of the array.

In this implementation, an array is a pointer to void.  The exact behavior of the
implementation is determined by structure TElementSpecification.  This, together
with the dimension of the array are stored before the first element of the array.
The dimension of the array can be obtained with the macro dw_DimA().

Additionally, mechanisms are provided for copying, initializing, printing, and
writing arrays.  If the array uses non-standard mechanisms for creating or 
destroying its elements, the functions dw_CopyArray() and dw_InitializeArray() 
should not be used.

=================================================================================
Specification of functions avaiable to arrays

Destructor
  void Destructor(void*)
    The destructor is called to destroy each element of the array.  It is assumed
    that each element of the array is in a valid state.  If the flag bit 
    dw_ARRAY_POINTER set, the calling syntax from an array is 

                          Destructor(((void**)a)[i])

    and the destructor should free the passed pointer if it is not null.  
    Otherwise the calling syntax is

                     Destructor((void*)(((char*)a) + i*size))

    and the destructor must not free the passed pointer.

---------------------------------------------------------------------------------

Default Constructor
  void DefaultConstructor(void *)
    The default constructor is called to initialize each element of the array.  
    The default constructor must not fail in the since that after a call to the 
    default constructor, the element is in a valid state.  The calling syntax is

                DefaultConstructor((void*)(((char*)a) + i*size))

    Note that calling syntax is the same for static and pointer arrays.  This 
    allows memory to be allocated to pointers.

---------------------------------------------------------------------------------

Pointer Copy Constructor
  void* PointerCopyConstructor(void*, void*)
    The pointer copy constructor is called when coping or initializing elements.  
    The calling syntax is 

                CopyConstructor(((void**)d)[i],((void**)s)[i])

    The pointer copy constructor returns a pointer to the destination.  If the
    destination was null, an attempt to create it can be made, and if successful
    a pointer to the newly allocated destination is returned.  Upon failure, a
    null pointer is returned.  In the event of failure, the pointer copy 
    constructor should leave the contents of the destination in a valid state.  
    It is assumed that source is in a valid state.

---------------------------------------------------------------------------------

Static Copy Constructor
  int StaticCopyConstructor(void*, void*)
    The static copy constructor is called when coping or initializing elements.  
    The calling syntax is 

     CopyConstructor((void*)(((char*)d) + i*size),(void*)(((char*)s) + i*size)

    The static copy constrctor should return one upon success and zero upon  
    failure.  In the event of failure, the static copy constructor should leave
    the contents of the destination in a state such that a call to the destructor
    will behave as expected.  The source is assume to be properly initialized.

---------------------------------------------------------------------------------
Print Routine
  int Print(FILE*, void*, char*)
    The print routine prints an element to the file f using the formating 
    information in the character string.  If the flag bit dw_ARRAY_POINTER is 
    set, the calling syntax is

                         Print(f,((void**)a)[i],format) 

    and otherwise is

                    Print(f,(void*)(((char*)a) + i*size),format)

---------------------------------------------------------------------------------

Read Routine
  int Read(FILE*, void*)
    The read routine reads an element from the file f.  If the flag bit
    dw_ARRAY_POINTER is set, the calling syntax is

                           Read(f,((void**)a)[i]) 

    and otherwise is

                      Read(f,(void*)(((char*)a) + i*size))

********************************************************************************/

#include <stdio.h>

//=========================== TElementSpecification ===========================//
#define dw_ARRAY_USE_MEMCPY      0x00000001
#define dw_ARRAY_POINTER         0x00000002
#define dw_ARRAY_ARRAY           0x00000004
#define dw_ARRAY_DELETE_SPECS    0x00000008

typedef struct
{
  int flag;

  int size;
  int offset;

  void (*destructor)(void*);
  void (*default_constructor)(void *);
  void* (*pointer_copy_constructor)(void*, void*);
  int (*static_copy_constructor)(void*, void*);
  int (*print_routine)(FILE*, void*, char*);
  int (*read_routine)(FILE*, void*);

} TElementSpecification;

TElementSpecification* CreateArraySpecification_pointer(void (*destructor)(void *));

extern TElementSpecification dw_IntSpecs;
extern TElementSpecification dw_DoubleSpecs;
extern TElementSpecification dw_FloatSpecs;
extern TElementSpecification dw_CharSpecs;
extern TElementSpecification dw_StringSpecs;
extern TElementSpecification dw_ArraySpecs;
extern TElementSpecification dw_PointerSpecs;
//=============================================================================//

//=== Macros ===
#define dw_DimA(a) (((int*)(a))[-1])
#define dw_SpecsA(a) (*((TElementSpecification**)(((char*)(a))-(sizeof(TElementSpecification*)+sizeof(int)))))
#define dw_IsArrayA(a) (dw_SpecsA(a)->flag & dw_ARRAY_ARRAY)

//=== Destructor ===//
void  dw_FreeArray(void* a);

//=== Constructors ===//
void* dw_CreateArray(TElementSpecification *specs, int dim);
void* dw_CreateMultidimensionalArray(TElementSpecification *specs, int depth, int *dim);
void* dw_CreateMultidimensionalArrayList(TElementSpecification *specs, int depth, ...);

//=== Routines ===//
void* dw_CopyArray(void* d, void* s);
int dw_PrintArray(FILE* f, void* a, char* format);
int dw_ReadArray(FILE* f, void* a);

// Array arrays
#define dw_CreateArray_array(dim) dw_CreateArray(&dw_ArraySpecs,dim)

// Pointer arrays
#define dw_CreateArray_pointer(dim,destructor)  dw_CreateArray(CreateArraySpecification_pointer(destructor),dim)
void DefaultPointerConstructor(void*);

// String arrays
#define dw_CreateArray_string(dim)  (char**)dw_CreateArray(&dw_StringSpecs,dim)
#define dw_CreateMultidimensionalArray_string(depth,dim) dw_CreateMultidimensionalArray(&dw_StringSpecs,depth,dim) 
void* dw_CreateMultidimensionalArrayList_string(int depth, ...);
#define dw_CreateRectangularArray_string(row,col) (char***)dw_CreateMultidimensionalArrayList_string(2,row,col)
#define dw_InitializeArray_string(a,x)  dw_InitializeArray(a,x)
               
// Integer arrays
#define dw_CreateArray_int(dim) (int*)dw_CreateArray(&dw_IntSpecs,dim)
#define dw_CreateMultidimensionalArray_int(depth,dim)  dw_CreateMultidimensionalArray(&dw_IntSpecs,depth,dim)
void* dw_CreateMultidimensionalArrayList_int(int depth, ...);
#define dw_CreateRectangularArray_int(row,col) (int**)dw_CreateMultidimensionalArrayList_int(2,row,col)
int dw_InitializeArray_int(void *a, int x);

// Double arrays
#define dw_CreateArray_double(dim) (double*)dw_CreateArray(&dw_DoubleSpecs,dim)
#define dw_CreateMultidimensionalArray_double(depth,dim) dw_CreateMultidimensionalArray(&dw_DoubleSpecs,depth,dim) 
void* dw_CreateMultidimensionalArrayList_double(int depth, ...);
#define dw_CreateRectangularArray_double(row,col) (double**)dw_CreateMultidimensionalArrayList_double(2,row,col)
int dw_InitializeArray_double(void *a, double x);

// Float arrays
#define dw_CreateArray_float(dim) (float*)dw_CreateArray(&dw_FloatSpecs,dim)
#define dw_CreateMultidimensionalArray_float(depth,dim) dw_CreateMultidimensionalArray(&dw_FloatSpecs,depth,dim)
void* dw_CreateMultidimensionalArrayList_float(int depth, ...);
#define dw_CreateRectangularArray_float(row,col) (float**)dw_CreateMultidimensionalArrayList_float(2,row,col)
int dw_InitializeArray_float(void *a, float x);

// Character arrays
#define dw_CreateArray_char(dim) (float*)dw_CreateArray(&dw_CharSpecs,dim)
#define dw_CreateMultidimensionalArray_char(depth,dim) dw_CreateMultidimensionalArray(&dw_CharSpecs,depth,dim)
void* dw_CreateMultidimensionalArrayList_char(int depth, ...);
#define dw_CreateRectangularArray_char(row,col) (char**)dw_CreateMultidimensionalArrayList_char(2,row,col)
int dw_InitializeArray_char(void *a, char x);

#endif
