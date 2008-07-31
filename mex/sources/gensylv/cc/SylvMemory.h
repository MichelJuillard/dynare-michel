/*
 * Copyright (C) 2003-2005 Ondra Kamenik
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvMemory.h,v 1.1.1.1 2004/06/04 13:00:49 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef SYLV_MEMORY_H
#define SYLV_MEMORY_H

#include "SylvParams.h"

#include <new>

class MallocAllocator {
#ifdef USE_MEMORY_POOL
public:
	void* operator new(size_t size);
	void* operator new[](size_t size);
	void operator delete(void* p);
	void operator delete[](void* p);
#endif
};

#ifdef USE_MEMORY_POOL
void* operator new(size_t size);
void* operator new[](size_t size);
void operator delete(void* p);
void operator delete[](void* p);
#endif

class SylvMemoryPool {
	char* base;
	size_t length;
	size_t allocated;
	bool stack_mode;
	SylvMemoryPool(const SylvMemoryPool&);
	const SylvMemoryPool& operator=(const SylvMemoryPool&);
public:
	SylvMemoryPool();
	~SylvMemoryPool();
	void init(size_t size);
	void* allocate(size_t size);
	void free(void* p);
	void reset();
	void setStackMode(bool);
};

class SylvMemoryDriver {
	SylvMemoryDriver(const SylvMemoryDriver&);
	const SylvMemoryDriver& operator=(const SylvMemoryDriver&);
public:
	SylvMemoryDriver(int num_d, int m, int n, int order);
	SylvMemoryDriver(const SylvParams& pars, int num_d, int m, int n, int order);
	static void setStackMode(bool);
	~SylvMemoryDriver();
protected:
	void allocate(int num_d, int m, int n, int order);
};

#endif /* SYLV_MEMORY_H */


// Local Variables:
// mode:C++
// End:
