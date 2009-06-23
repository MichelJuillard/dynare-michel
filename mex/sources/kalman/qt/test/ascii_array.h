// ascii_matrix.h 
// Based on work of  Ondra Kamenik


struct AsciiNumberArray {
	int rows;
	int cols;
  double * data;
	char* fname;
  void GetMX(const char* fname,  int rows, int cols);
  void WriteMX();
};

void WriteMX(char* fname, double* data, int rows, int cols);

