/* $Header: /var/lib/cvs/dynare_cpp/sylv/testing/MMMatrix.cpp,v 1.1.1.1 2004/06/04 13:01:13 kamenik Exp $ */

/* Tag $Name:  $ */

#include "MMMatrix.h"

#include <cstdio>
#include <cstring>

MMMatrixIn::MMMatrixIn(const char* fname)
{
	FILE* fd;
	if (NULL==(fd = fopen(fname,"r")))
		throw MMException(string("Cannot open file ")+fname+" for reading\n");

	char buffer[1000];
	// jump over initial comments
	while (fgets(buffer, 1000, fd) && strncmp(buffer, "%%", 2)) {}
	// read in number of rows and cols
	if (!fgets(buffer, 1000, fd))
		throw MMException(string("Cannot read rows and cols while reading ")+fname+"\n");
	if (2 != sscanf(buffer, "%d %d", &rows, &cols))
		throw MMException("Couldn't parse rows and cols\n");
	// read in data
	data = (double*) operator new[](rows*cols*sizeof(double));
	int len = rows*cols;
	int i = 0;
	while (fgets(buffer, 1000, fd) && i < len) {
		if (1 != sscanf(buffer, "%lf", &data[i]))
			throw MMException(string("Couldn't parse float number ")+buffer+"\n");
		i++;
	}
	if (i < len) {
		char mes[1000];
		sprintf(mes,"Couldn't read all %d lines, read %d so far\n",len,i);
		throw MMException(mes);
	}
	fclose(fd);
}

MMMatrixIn::~MMMatrixIn()
{
	operator delete [](data);
}


void MMMatrixOut::write(const char* fname, int rows, int cols, const double* data)
{
	FILE* fd;
	if (NULL==(fd = fopen(fname,"w")))
		throw MMException(string("Cannot open file ")+fname+" for writing\n");

	if (0 > fprintf(fd, "%%%%MatrixMarket matrix array real general\n"))
		throw MMException(string("Output error when writing file ")+fname);
	if (0 > fprintf(fd, "%d %d\n", rows, cols))
		throw MMException(string("Output error when writing file ")+fname);
	int running = 0;
	for (int i = 0; i < cols; i++) {
		for (int j = 0 ; j < rows; j++) {
			if (0 > fprintf(fd, "%40.35g\n", data[running]))
				throw MMException(string("Output error when writing file ")+fname);
			running++;
		}
	}
	fclose(fd);
}

void MMMatrixOut::write(const char* fname, const GeneralMatrix& m)
{
	write(fname, m.numRows(), m.numCols(), m.base());
}
