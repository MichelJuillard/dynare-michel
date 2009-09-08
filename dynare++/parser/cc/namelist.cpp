// Copyright (C) 2006, Ondra Kamenik

// $Id: namelist.cpp 42 2007-01-22 21:53:24Z ondra $

#include "namelist.h"

#include <string.h>

using namespace ogp;

/** A global symbol for passing info to NameListParser from its
 * parser. */
NameListParser* name_list_parser;

void* namelist__scan_buffer(char*, unsigned int);
void namelist__destroy_buffer(void*);
void namelist_parse();

void NameListParser::namelist_parse(int length, const char* stream)
{
	char* buffer = new char[length+2];
	strncpy(buffer, stream, length);
	buffer[length] = '\0';
	buffer[length+1] = '\0';
	void* p = namelist__scan_buffer(buffer, (unsigned int)length+2);
	name_list_parser = this;
	::namelist_parse();
	delete [] buffer;
	namelist__destroy_buffer(p);
}
