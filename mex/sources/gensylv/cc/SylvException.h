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

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvException.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef SYLV_EXCEPTION_H
#define SYLV_EXCEPTION_H

#include "SylvMemory.h"


class SylvException : public MallocAllocator {
protected:
	char file[50];
	int line;
	const SylvException* source;
public:
	SylvException(const char* f, int l, const SylvException* s);
	virtual ~SylvException();
	virtual int printMessage(char* str, int maxlen) const;
	void printMessage() const;
};

class SylvExceptionMessage : public SylvException {
	char message[500];
public:
	SylvExceptionMessage(const char* f, int l, const char* mes);
	virtual int printMessage(char* str, int maxlen) const;
};

// define macros:
#define SYLV_EXCEPTION(exc) (SylvException(__FILE__, __LINE__, exc))
#define SYLV_MES_EXCEPTION(mes) (SylvExceptionMessage(__FILE__, __LINE__, mes)) 

#endif /* SYLV_EXCEPTION_H */


// Local Variables:
// mode:C++
// End:
