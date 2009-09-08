// Copyright (C) 2006, Ondra Kamenik

// $Id: dynare_exception.h 853 2006-08-01 08:42:42Z kamenik $

#ifndef DYNARE_EXCEPTION_H
#define DYNARE_EXCEPTION_H

#include <string>

class DynareException {
	char* mes;
public:
	DynareException(const char* m, const char* fname, int line, int col)
		{
			mes = new char[strlen(m) + strlen(fname) + 100];
			sprintf(mes, "Parse error at %s, line %d, column %d: %s", fname, line, col, m);
		}
	DynareException(const char* fname, int line, const std::string& m)
		{
			mes = new char[m.size() + strlen(fname) + 50];
			sprintf(mes, "%s:%d: %s", fname, line, m.c_str());
		}
	DynareException(const char* m, int offset)
		{
			mes = new char[strlen(m) + 100];
			sprintf(mes, "Parse error in provided string at offset %d: %s", offset, m);
		}
	DynareException(const DynareException& e)
		: mes(new char[strlen(e.mes)+1])
		{strcpy(mes, e.mes);}
	virtual ~DynareException()
		{delete [] mes;}
	const char* message() const
		{return mes;}
};

#endif

// Local Variables:
// mode:C++
// End:
