/*1:*/

#ifndef TL_EXCEPTION_H
#define TL_EXCEPTION_H

#include <string.h> 
#include <stdio.h> 

/*2:*/

#ifndef TL_DEBUG_EXCEPTION
#define TL_DEBUG_EXCEPTION 1
#endif

#ifndef TL_DEBUG
#define TL_DEBUG 0
#endif

#define TL_RAISE(mes) \
if (TL_DEBUG >= TL_DEBUG_EXCEPTION) throw TLException(__FILE__, __LINE__, mes);

#define TL_RAISE_IF(expr, mes) \
if (TL_DEBUG >= TL_DEBUG_EXCEPTION && (expr)) throw TLException(__FILE__, __LINE__, mes);

/*3:*/

class TLException{
	char fname[50];
	int lnum;
	char message[500];
public:
	TLException(const char*f,int l,const char*mes)
	{
		strncpy(fname,f,50);fname[49]= '\0';
		strncpy(message,mes,500);message[499]= '\0';
		lnum= l;
	}
	virtual~TLException(){}
	virtual void print()const
	{printf("At %s:%d:%s\n",fname,lnum,message);}
};


/*:3*/
;

/*:2*/
;

#endif

/*:1*/
