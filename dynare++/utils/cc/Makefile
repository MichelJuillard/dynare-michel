# Copyright (C) 2005, Ondra Kamenik

# $Id: Makefile 843 2006-07-28 08:54:19Z tamas $

CC_FLAGS := -Wall -I../../

ifeq ($(DEBUG),yes)
	CC_FLAGS := $(CC_FLAGS) -g
else
	CC_FLAGS := $(CC_FLAGS) -O3
endif

objects := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
headers := $(wildcard *.h)
source  := $(wildcard *.cpp)

%.o: %.cpp $(headers)
	$(CC) $(CC_FLAGS) -c $*.cpp

utils.a: $(objects) $(headers) $(source)
	ar cr utils.a $(objects)

clear:
	rm -f *~
	rm -f *.o
	rm -f utils.a
