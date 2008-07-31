/*
 * Copyright (C) 2003-2005 Ondra Kamenik
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SylvException.cpp,v 1.2 2004/10/01 10:30:40 kamenik Exp $ */

/* Tag $Name:  $ */

#include "SylvException.h"

#include <string.h>
#include <stdio.h>

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
	printf(mes);
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


