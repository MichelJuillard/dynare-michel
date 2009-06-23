// ascii_matrix.cpp 
// Based on work of 2005, Ondra Kamenik

#include "ascii_array.h"

#include <stdio.h>
#include <string.h>
#include <iostream>

#include <fstream>
#include <string>

// if the file doesn't exist, the number array is empty
void 
AsciiNumberArray::GetMX(const char* fname, int INrows, int INcols)
  {
  rows = 0;
  cols = INcols;
  
  std::ifstream file(fname);
  std::string line;
  
  data=(double*)calloc(INcols*INrows,sizeof(double));
  while (getline(file, line)) 
    {
    rows++;
    int icols = 0;
    const char delims[] = " \r\n\t";
    char* lineptr = strdup(line.c_str());
    char* tok = strtok(lineptr, delims);
    while (tok) 
      {
      icols++;
      double item;
      if (1 != sscanf(tok, "%lf", &item)) 
        {
        fprintf(stderr, "Couldn't parse a token %s as double.\n", tok);
        exit(1);
        }
      data[(rows-1)*INcols+icols-1]=item;
      tok = strtok(NULL, delims);
      }
    free(lineptr);
    if (cols) 
      {
      if (cols != icols) 
        {
        fprintf(stderr, "Asserted a different number of columns.\n");
        exit(1);
        }
      } 
    else 
      {
      cols = icols;
      }
    }
  }


void 
AsciiNumberArray::WriteMX()
  {
  std::ofstream outFile(strcat(fname, "_out"));
  for (int i = 0; i < rows; i++)
    {
    for (int j = 0; j < cols; j++)
      {
      outFile << data[i*cols+j] << "\t";
      }
    outFile << std::endl;
    }
  outFile.close();
  }


void WriteMX(char* fname, double* data, int rows, int cols)
  {
  AsciiNumberArray OutArray;
  OutArray.fname=fname;
  OutArray.rows=rows;
  OutArray.cols=cols;
  OutArray.data=data;
  OutArray.WriteMX();
  }
