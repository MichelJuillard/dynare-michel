// $Id: ascii_matrix.h 534 2005-11-30 13:58:11Z kamenik $
// Copyright 2005, Ondra Kamenik

#include "GeneralMatrix.h"

#include <vector>
#include <string>

struct AsciiNumberArray {
	int rows;
	int cols;
	std::vector<double> data;
	AsciiNumberArray(const char* fname)
		{parse(fname);}
	AsciiNumberArray(std::string fname)
		{parse(fname.c_str());}
protected:
	void parse(const char* fname);
};

class AsciiMatrix : public GeneralMatrix {
public:
	AsciiMatrix(const AsciiNumberArray& na);
};
