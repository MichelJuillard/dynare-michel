// $Id: ascii_matrix.cpp 534 2005-11-30 13:58:11Z kamenik $
// Copyright 2005, Ondra Kamenik

#include "ascii_matrix.h"

#include <stdio.h>
#include <string.h>

#include <fstream>
#include <string>

// if the file doesn't exist, the number array is empty
void AsciiNumberArray::parse(const char* fname)
{
	rows = 0;
	cols = 0;

	std::ifstream file(fname);
	std::string line;

	while (getline(file, line)) {
		rows++;
		int icols = 0;
		const char delims[] = " \r\n\t";
		char* lineptr = strdup(line.c_str());
		char* tok = strtok(lineptr, delims);
		while (tok) {
			icols++;
			double item;
			if (1 != sscanf(tok, "%lf", &item)) {
				fprintf(stderr, "Couldn't parse a token %s as double.\n", tok);
				exit(1);
			}
			data.push_back(item);
			tok = strtok(NULL, delims);
		}
		free(lineptr);
		if (cols) {
			if (cols != icols) {
				fprintf(stderr, "Asserted a different number of columns.\n");
				exit(1);
			}
		} else {
			cols = icols;
		}
	}
}


AsciiMatrix::AsciiMatrix(const AsciiNumberArray& na)
	: GeneralMatrix(na.rows, na.cols)
{
	for (int i = 0; i < numRows(); i++)
		for (int j = 0; j < numCols(); j++)
			get(i, j) = na.data[i*numCols()+j];
}
