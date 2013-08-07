// Copyright (C) 2006, Ondra Kamenik

// $Id: matrix_parser.cpp 2269 2008-11-23 14:33:22Z michel $

#include "parser_exception.h"
#include "matrix_parser.h"
#include "location.h"
#include "matrix_tab.hh"
#include <cstring>

using namespace ogp;

/** A global symbol for passing info to the MatrixParser from
 * matrix_parse(). */
MatrixParser* mparser;

/** The declaration of functions defined in matrix_ll.cc and
 * matrix_tab.cc generated from matrix.lex and matrix.y. */
void* matrix__scan_buffer(char*, size_t);
void matrix__destroy_buffer(void*);
int matrix_parse();
extern ogp::location_type matrix_lloc;

void MatrixParser::parse(int length, const char* stream)
{
	// reinitialize the object
	data.clear();
	row_lengths.clear();
	nc = 0;
	// allocate temporary buffer and parse
	char* buffer = new char[length+2];
	strncpy(buffer, stream, length);
	buffer[length] = '\0';
	buffer[length+1] = '\0';
	matrix_lloc.off = 0;
	matrix_lloc.ll = 0;
	void* p = matrix__scan_buffer(buffer, (unsigned int)length+2);
	mparser = this;
	matrix_parse();
	delete [] buffer;
	matrix__destroy_buffer(p);
}

void MatrixParser::add_item(double v)
{
	data.push_back(v);
	if (row_lengths.size() == 0)
		row_lengths.push_back(0);
	(row_lengths.back())++;
	if (row_lengths.back() > nc)
		nc = row_lengths.back();
}

void MatrixParser::start_row()
{
	row_lengths.push_back(0);
}

void MatrixParser::error(const char* mes) const
{
	throw ParserException(mes, matrix_lloc.off);
}

int MatrixParser::find_first_non_empty_row(int start) const
{
	int r = start;
	while (r < (int)row_lengths.size() && row_lengths[r] == 0)
		r++;
	return r;
}

MPIterator MatrixParser::begin() const
{
	MPIterator it(*this);
	return it;
}

MPIterator MatrixParser::end() const
{
	MPIterator it(*this, "end");
	return it;
}

MPIterator::MPIterator(const MatrixParser& mp)
	: p(&mp), i(0), c(0), r(mp.find_first_non_empty_row())
{}

MPIterator::MPIterator(const MatrixParser& mp, const char* dummy)
	: p(&mp), i(mp.data.size()), c(0), r(mp.row_lengths.size())
{}

MPIterator& MPIterator::operator++()
{
	i++;
	c++;
	if (p->row_lengths[r] <= c) {
		c = 0;
		r = p->find_first_non_empty_row(r+1);
	}
	return *this;
}
