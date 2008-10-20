/*1:*/

#ifndef KORD_EXCEPTION_H
#define KORD_EXCEPTION_H

#include <string.h> 
#include <stdio.h> 

#define KORD_RAISE(mes) \
throw KordException(__FILE__, __LINE__, mes);

#define KORD_RAISE_IF(expr, mes) \
if (expr) throw KordException(__FILE__, __LINE__, mes);

#define KORD_RAISE_X(mes, c) \
throw KordException(__FILE__, __LINE__, mes, c);

#define KORD_RAISE_IF_X(expr, mes, c) \
if (expr) throw KordException(__FILE__, __LINE__, mes, c);

/*2:*/

class KordException{
protected:
	char fname[50];
	int lnum;
	char message[500];
	int cd;
public:
	KordException(const char*f,int l,const char*mes,int c= 255)
	{
		strncpy(fname,f,50);fname[49]= '\0';
		strncpy(message,mes,500);message[499]= '\0';
		lnum= l;
		cd= c;
	}
	virtual~KordException(){}
	virtual void print()const
	{printf("At %s:%d:(%d):%s\n",fname,lnum,cd,message);}
	virtual int code()const
	{return cd;}
	const char*get_message()const
	{return message;}
};

/*:2*/
;
/*3:*/

#define KORD_FP_NOT_CONV 254
#define KORD_FP_NOT_FINITE 253
#define KORD_MD_NOT_STABLE 252

/*:3*/
;

#endif

/*:1*/
