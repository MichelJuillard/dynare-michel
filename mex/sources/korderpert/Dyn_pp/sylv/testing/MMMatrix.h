/* $Header: /var/lib/cvs/dynare_cpp/sylv/testing/MMMatrix.h,v 1.1.1.1 2004/06/04 13:01:13 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef MM_MATRIX_H
#define MM_MATRIX_H

#include "GeneralMatrix.h"
#include "SylvMemory.h"

#include <string>

using namespace std;

class MMException : public MallocAllocator {
	string message;
public:
	MMException(string mes) : message(mes) {}
	MMException(const char* mes) : message(mes) {}
	const char* getMessage() const {return message.data();}
};

class MMMatrixIn : public MallocAllocator {
	double* data;
	int rows;
	int cols;
public:
	MMMatrixIn(const char* fname);
	~MMMatrixIn();
	const double* getData() const {return data;}
	int size() const {return rows*cols;}
	int row() const {return rows;}
	int col() const {return cols;}
};

class MMMatrixOut : public MallocAllocator {
public:
	static void write(const char* fname, int rows, int cols, const double* data);
	static void write(const char* fname, const GeneralMatrix& m);
};

#endif /* MM_MATRIX_H */

// Local Variables:
// mode:C++
// End:
