/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvException.cpp,v 1.2 2004/10/01 10:30:40 kamenik Exp $ */

/* Tag $Name:  $ */

#include "SylvException.h"

#include <cstring>
#include <cstdio>

SylvException::SylvException(const char* f, int l, const SylvException* s)
{
	strcpy(file,f);
	line = l;
	source = s;
}

SylvException::~SylvException()
{
	if (source != NULL) {
		delete source;
	}
}

void SylvException::printMessage() const
{
	char mes[1500];
	mes[0] = '\0';
	printMessage(mes, 1499);
	puts(mes);
}

int SylvException::printMessage(char* str, int maxlen) const
{
	int remain = maxlen;
	if (source != NULL) {
		remain = source->printMessage(str, maxlen);
	}
	char aux[100];
	sprintf(aux, "From %s:%d\n", file, line);
	int newremain = remain - strlen(aux);
	if (newremain < 0) {
		aux[remain] = '\0';
		newremain = 0;
	}
	strcat(str, aux);
	return newremain;
}

SylvExceptionMessage::SylvExceptionMessage(const char* f, int i,
										   const char* mes)
	: SylvException(f,i,NULL)
{
	strcpy(message,mes);
}

int SylvExceptionMessage::printMessage(char* str, int maxlen) const
{
	char aux[600];
	sprintf(aux, "At %s:%d:%s\n", file, line, message);
	int newremain = maxlen - strlen(aux);
	if (newremain < 0) {
		aux[maxlen] = '\0';
		newremain = 0;
	}
	strcat(str, aux);
	return newremain;
}


