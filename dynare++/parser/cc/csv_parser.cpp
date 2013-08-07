#include "csv_parser.h"
#include "parser_exception.h"
#include "location.h"
#include "csv_tab.hh"
#include <cstring>

using namespace ogp;

/** A global symbol for passing info to the CSVParser from
 * csv_parse(). */
CSVParser* csv_parser;

/** The declaration of functions defined in csv_ll.cc and
 * csv_tab.cc generated from csv.lex and csv.y. */
void* csv__scan_buffer(char*, unsigned int);
void csv__destroy_buffer(void*);
int csv_parse();

extern ogp::location_type csv_lloc;

void CSVParser::csv_error(const char* mes)
{
	throw ParserException(mes, csv_lloc.off);
}

void CSVParser::csv_parse(int length, const char* str)
{
	// allocate temporary buffer and parse
	char* buffer = new char[length+2];
	strncpy(buffer, str, length);
	buffer[length] = '\0';
	buffer[length+1] = '\0';
	csv_lloc.off = 0;
	csv_lloc.ll = 0;
	parsed_string = buffer;
	void* p = csv__scan_buffer(buffer, (unsigned int)length+2);
	csv_parser = this;
	::csv_parse();
	delete [] buffer;
	csv__destroy_buffer(p);
	parsed_string = NULL;
}
