// Copyright (C) 2005, Ondra Kamenik

// $Id: exception.h 1367 2007-07-11 14:21:57Z kamenik $

#ifndef OGU_EXCEPTION_H
#define OGU_EXCEPTION_H

#include <cstdio>
#include <cstring>

#include <string>
#include <algorithm>

namespace ogu {

	/** A primitive exception. */
	class Exception {
		static const int file_length = 100;
		static const int mes_length = 500;
	protected:
		char file[file_length];
		int line;
		char mes[mes_length];
	public:
		Exception(const char* f, int l, const char* m)
			{
				strncpy(file, f, file_length-1);
				file[file_length-1] = '\0';
				line = l;
#ifndef _MSC_VER
				strncpy(mes, m, std::min(mes_length-1,(int)strlen(m)));
#else
				/* MSVC doesn't define std::min() (should be in <algorithm>),
				   but instead has a macro in <windows.h> */
				strncpy(mes, m, min(mes_length-1,(int)strlen(m)));
#endif
				mes[mes_length-1] = '\0';
			}
		Exception(const char* f, int l, const std::string& m)
			{
				strncpy(file, f, file_length-1);
				file[file_length-1] = '\0';
				line = l;
#ifndef _MSC_VER
				strncpy(mes, m.c_str(), std::min(mes_length-1,(int)m.length()));
#else
				/* MSVC doesn't define std::min() (should be in <algorithm>),
				   but instead has a macro in <windows.h> */
				strncpy(mes, m.c_str(), min(mes_length-1,(int)m.length()));
#endif
				mes[mes_length-1] = '\0';
			}
		virtual ~Exception() {}
		void print(FILE* fd) const
			{ fprintf(fd, "%s:%d: %s\n", file, line, mes); }
		void print() const
			{ print(stdout); }
		const char* message() const
			{ return mes; }
	};
};

#endif

// Local Variables:
// mode:C++
// End:
