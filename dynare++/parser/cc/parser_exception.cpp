// Copyright (C) 2006, Ondra Kamenik

// $Id: parser_exception.cpp 2269 2008-11-23 14:33:22Z michel $

#include "parser_exception.h"
#include <cstring>

using namespace ogp;

ParserException::ParserException(const char* m, int offset)
	: mes(new char[strlen(m)+1]), off(offset),
	  aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
	strcpy(mes, m);
}

ParserException::ParserException(const string& m, int offset)
	: mes(new char[m.size()+1]), off(offset),
	  aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
	strncpy(mes, m.c_str(), m.size());
	mes[m.size()] = '\0';
}

ParserException::ParserException(const string& m, const char* dum, int i1)
	: mes(new char[m.size()+1]), off(0),
	  aux_i1(i1), aux_i2(-1), aux_i3(-1)
{
	strncpy(mes, m.c_str(), m.size());
	mes[m.size()] = '\0';
}

ParserException::ParserException(const string& m, const char* dum, int i1, int i2) 
	: mes(new char[m.size()+1]), off(0),
	  aux_i1(i1), aux_i2(i2), aux_i3(-1)
{
	strncpy(mes, m.c_str(), m.size());
	mes[m.size()] = '\0';
}

ParserException::ParserException(const string& m, const char* dum, int i1, int i2, int i3) 
	: mes(new char[m.size()+1]), off(0),
	  aux_i1(i1), aux_i2(i2), aux_i3(i3)
{
	strncpy(mes, m.c_str(), m.size());
	mes[m.size()] = '\0';
}

ParserException::ParserException(const ParserException& m, int plus_offset)
	: mes(NULL),
	  aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
	copy(m);
	off += plus_offset;
}

ParserException::ParserException(const ParserException& m, const char* dum, int i)
	: mes(NULL),
	  aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
	copy(m);
	aux_i3 = m.aux_i2;
	aux_i2 = m.aux_i1;
	aux_i1 = i;
}

ParserException::ParserException(const ParserException& m, const char* dum, int i1, int i2)
	: mes(NULL),
	  aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
	copy(m);
	aux_i3 = m.aux_i1;
	aux_i2 = i2;
	aux_i1 = i1;
}

ParserException::ParserException(const ParserException& m, const char* dum, int i1, int i2, int i3)
	: mes(NULL),
	  aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
	copy(m);
	aux_i3 = i3;
	aux_i2 = i2;
	aux_i1 = i1;
}


ParserException::ParserException(const ParserException& e)
	: mes(NULL),
	  aux_i1(-1), aux_i2(-1), aux_i3(-1)
{
	copy(e);
} 

ParserException::~ParserException()
{
	delete [] mes;
}

void ParserException::copy(const ParserException& e)
{
	if (mes)
		delete [] mes;
	mes = new char[strlen(e.mes)+1];
	strcpy(mes, e.mes);
	off = e.off;
	aux_i1 = e.aux_i1;
	aux_i2 = e.aux_i2;
	aux_i3 = e.aux_i3;
}

void ParserException::print(FILE* fd) const
{
	// todo: to be refined
	fprintf(fd, "%s: offset %d\n", mes, off);
}
