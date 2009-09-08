// Copyright (C) 2006, Ondra Kamenik

// $Id: location.h 762 2006-05-22 13:00:07Z kamenik $

// Purpose: This file defines macros for lex and bison so that the
// very primitive location tracking would be enabled. The location of
// a token is given by offset of its first character. The offset is
// relative to the number which is (and must be) initialized before
// parsing. This file is to be included to the top of bison and lex
// sources.

// How to use: in preamble of bison and flex, you must include this
// file and declare extern YYLTYPE prefix##lloc. In addition, in flex,
// you must define int prefix##ll =0; and use macro SET_LLOC(prefix)
// in EVERY action consuming material (this can be done with #define
// YY_USER_ACTION) and in bison you must use option %locations.


#ifndef OG_LOCATION_H
#define OG_LOCATION_H

namespace ogp {

	struct location_type {
		int off; // offset of the token
		int ll; // length ot the token
		location_type() : off(0), ll(0) {}
	};

};

#define YYLTYPE ogp::location_type

// set current off to the first off and add all lengths
#define YYLLOC_DEFAULT(Current, Rhs, N) \
  {(Current).off    =  (Rhs)[1].off;    \
   (Current).ll     =  0;               \
   for (int i = 1; i <= N; i++) (Current).ll += (Rhs)[i].ll;}

#define SET_LLOC(prefix) (prefix##lloc.off += prefix##lloc.ll, prefix##lloc.ll = prefix##leng)

#endif

// Local Variables:
// mode:C++
// End:
