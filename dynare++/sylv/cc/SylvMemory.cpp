/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvMemory.cpp,v 1.1.1.1 2004/06/04 13:00:49 kamenik Exp $ */

/* Tag $Name:  $ */

#include "SylvMemory.h"
#include "SylvException.h"
#include "KronVector.h"

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
# include <dynmex.h>
#endif

#include <math.h> 
#include <stdio.h>
#include <stdlib.h>

/**********************************************************/
/*   SylvMemoryPool                                       */
/**********************************************************/

SylvMemoryPool memory_pool;

SylvMemoryPool::SylvMemoryPool()
	: base(0), length(0), allocated(0), stack_mode(false)
{
}

void SylvMemoryPool::init(size_t size)
{
#ifdef USE_MEMORY_POOL
	length = size;

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
	if (base)
		throw SYLV_MES_EXCEPTION("Attempt to use matlab memory pool twice.");
	base = (char*) mxMalloc(length);
#else
	base = (char*) malloc(length);
#endif

#else
	throw SYLV_MES_EXCEPTION("SylvMemoryPool::init() called for non memory pool code.");
#endif
}

void* SylvMemoryPool::allocate(size_t size)
{
#ifdef USE_MEMORY_POOL
	if (allocated + size < length) {
		char* res = base + allocated;
		allocated += size;
		return res;
	} else {
		throw SYLV_MES_EXCEPTION("Run out of memory space");
	}
#else
	throw SYLV_MES_EXCEPTION("SylvMemoryPool::allocate() called for non memory pool code.");
#endif
}

void SylvMemoryPool::free(void* p)
{
#ifdef USE_MEMORY_POOL
	int offset = ((char*)p) - base;

#ifdef DEBUG
	if (offset < 0)
		throw SYLV_MES_EXCEPTION("SylvMemoryPool::free() frees wrong address < begin.");
	if (offset >= (int)length)
		throw SYLV_MES_EXCEPTION("SylvMemoryPool::free() frees wrong address > end.");
#endif	

	if (stack_mode && offset >= 0 && offset < (int)allocated)
		allocated = offset;

#else
	throw SYLV_MES_EXCEPTION("SylvMemoryPool::free() called for non memory pool code.");
#endif
}

void SylvMemoryPool::setStackMode(bool mode)
{
	stack_mode = mode;
}

SylvMemoryPool::~SylvMemoryPool()
{
	reset();
}

void SylvMemoryPool::reset()
{
#if !defined(MATLAB_MEX_FILE) && !defined(OCTAVE_MEX_FILE)
	delete [] base;
	base = 0;
	allocated = 0;
	length = 0;
	stack_mode = false;
#endif
}

/**********************************************************/
/*   global new and delete                                */
/**********************************************************/

#ifdef USE_MEMORY_POOL

void* operator new(size_t size)
{
	return memory_pool.allocate(size);
}

void* operator new[](size_t size)
{
	return memory_pool.allocate(size);
}

void operator delete(void* p)
{
	memory_pool.free(p);
}

void operator delete[](void* p)
{
	memory_pool.free(p);
}

#endif

/**********************************************************/
/*   saved version of global new and delete               */
/**********************************************************/

#ifdef USE_MEMORY_POOL
void* MallocAllocator::operator new(size_t size)
{
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
	throw SYLV_MES_EXCEPTION("Attempt to call wrong memory allocator.");
#else
	void* res = malloc(size);
	if (!res) 
		throw SYLV_MES_EXCEPTION("Malloc unable to allocate memory.");
	return res;
#endif
}

void* MallocAllocator::operator new[](size_t size)
{
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
	throw SYLV_MES_EXCEPTION("Attempt to call wrong memory allocator.");
#else
	void* res = malloc(size);
	if (!res)
		throw SYLV_MES_EXCEPTION("Malloc unable allocate memory.");
	return res;
#endif
}

void MallocAllocator::operator delete(void* p)
{
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
	throw SYLV_MES_EXCEPTION("Attempt to call wrong memory destructor.");
#else
	free(p);
#endif
}

void MallocAllocator::operator delete[](void* p)
{
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
	throw SYLV_MES_EXCEPTION("Attempt to call wrong memory destructor.");
#else
	free(p);
#endif
}

#endif


/**********************************************************/
/*   SylvMemoryDriver                                     */
/**********************************************************/

void SylvMemoryDriver::allocate(int num_d, int m, int n, int order)
{
#ifdef USE_MEMORY_POOL
	int x_cols = power(m,order);
	int total = num_d*x_cols*n; // storage for big matrices
	total += x_cols; // storage for one extra row of a big matrix
	int dig_vectors = (int)ceil(((double)(power(m,order)-1))/(m-1));
	total += 8*n*dig_vectors; // storage for kron vectors instantiated during solv
	total += 50*(m*m+n*n); // some storage for small square matrices
	total *= sizeof(double); // everything in doubles
	memory_pool.init(total);
#endif
}


SylvMemoryDriver::SylvMemoryDriver(int num_d, int m, int n, int order)
{
	allocate(num_d, m, n, order);
}

SylvMemoryDriver::SylvMemoryDriver(const SylvParams& pars, int num_d,
								   int m, int n, int order)
{
	if (*(pars.method) == SylvParams::iter)
		num_d++;
	if (*(pars.want_check))
		num_d++;
	allocate(num_d, m, n, order);
}

SylvMemoryDriver::~SylvMemoryDriver()
{
	memory_pool.reset();
}

void SylvMemoryDriver::setStackMode(bool mode) {
	memory_pool.setStackMode(mode);
}
