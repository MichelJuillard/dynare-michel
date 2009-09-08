// Copyright (C) 2006, Ondra Kamenik

// $Id: parser_exception.h 1761 2008-03-31 14:27:13Z kamenik $

#ifndef OG_FORMULA_PARSER_H
#define OG_FORMULA_PARSER_H

#include <string>

namespace ogp {
	using std::string;

	/** This is an easy exception, which, besides the message, stores
	 * also an offset of the parse error. Since we might need to track
	 * the argument number and for example the filed in the argument
	 * which caused the error, we add three integers, which have no
	 * semantics here. They should be documented in the function which
	 * throws an exception and sets them. Their default value is -1,
	 * which means they have not been set. */
	class ParserException {
	protected:
		char* mes;
		int off;
		int aux_i1;
		int aux_i2;
		int aux_i3;
	public:
		ParserException(const char* m, int offset);
		ParserException(const string& m, int offset);
		ParserException(const string& m, const char* dum, int i1);
		ParserException(const string& m, const char* dum, int i1, int i2);
		ParserException(const string& m, const char* dum, int i1, int i2, int i3);
		ParserException(const ParserException& e, int plus_offset);
		/** Makes a copy and pushes given integer to aux_i1 shuffling
		 * others and forgetting the last. */
		ParserException(const ParserException& e, const char* dum, int i);
		/** Makes a copy and pushes given two integers to aux_i1 and aux_i2  shuffling
		 * others and forgetting the last two. */
		ParserException(const ParserException& e, const char* dum, int i1, int i2);
		/** Makes a copy and pushes given three integers to aux_i1, aux_i2, aus_i3 shuffling
		 * others and forgetting the last three. */
		ParserException(const ParserException& e, const char* dum, int i1, int i2, int i3);
		ParserException(const ParserException& e);
		virtual ~ParserException();
		void print(FILE* fd) const;
		const char* message() const
			{return mes;}
		int offset() const
			{return off;}
		const int& i1() const
			{return aux_i1;}
		int& i1()
			{return aux_i1;}
		const int& i2() const
			{return aux_i2;}
		int& i2()
			{return aux_i2;}
		const int& i3() const
			{return aux_i3;}
		int& i3()
			{return aux_i3;}
	protected:
		void copy(const ParserException& e);
	};
};

#endif

// Local Variables:
// mode:C++
// End:
