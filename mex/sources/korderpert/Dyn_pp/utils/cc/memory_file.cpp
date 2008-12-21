// Copyright (C) 2005, Ondra Kamenik

// $Id: memory_file.cpp 987 2006-10-17 14:39:19Z kamenik $

#include "memory_file.h"

#include <stdio.h>

using namespace ogu;

int ogu::calc_pos_offset(int length, const char* str, int line, int col)
{
	int i = 0;
	int il = 1;
	int ic = 1;
	while (i < length && il <= line && ic <= col) {
		if (str[i] == '\n') {
			il++;
			ic = 1;
		} else {
			ic++;
		}
	}
	return i;
}

void ogu::calc_pos_line_and_col(int length, const char* str, int offset,
						   int& line, int& col)
{
	line = 1;
	col = 0;
	int i = 0;
	while (i < length && i < offset) {
		if (str[i] == '\n') {
			line++;
			col = 0;
		}
		i++;
		col++;
	}
}

MemoryFile::MemoryFile(const char* fname)
	: len(-1), data(NULL)
{
	FILE* fd = fopen(fname, "rb");
	if (fd) {
		// get the file size
		fseek(fd, 0, SEEK_END);
		len = ftell(fd);
		// allocate space for the file plus ending '\0' character
		data = new char[len+1];
		// read file and set data
		fseek(fd, 0, SEEK_SET);
		int i = 0;
		int c;
		while (EOF != (c = fgetc(fd)))
			data[i++] = (unsigned char)c;
		data[len] = '\0';
		fclose(fd);
	}
}
