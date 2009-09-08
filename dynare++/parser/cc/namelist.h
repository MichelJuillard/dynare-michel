// Copyright (C) 2007, Ondra Kamenik

// $Id: namelist.h 107 2007-05-10 22:35:04Z ondra $

#ifndef OGP_NAMELIST
#define OGP_NAMELIST

namespace ogp {

	/** Parent class of all parsers parsing a namelist. They must
	 * implement add_name() method and error() method, which is called
	 * when an parse error occurs. 
	 *
	 * Parsing a name list is done as follows: implement
	 * NameListParser interface, create the object, and call
	 * NameListParser::namelist_parse(int lengt, const char*
	 * text). When implementing error(), one may consult global
	 * location_type namelist_lloc. */
	class NameListParser {
	public:
		virtual ~NameListParser() {}
		virtual void add_name(const char* name) = 0;
		virtual void namelist_error(const char* mes) = 0;
		void namelist_parse(int length, const char* text);
	};
};

#endif

// Local Variables:
// mode:C++
// End:
