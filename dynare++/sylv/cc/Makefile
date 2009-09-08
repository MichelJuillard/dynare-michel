# $Header: /var/lib/cvs/dynare_cpp/sylv/cc/Makefile,v 1.4 2005/01/18 21:28:26 kamenik Exp $

# Tag $Name:  $

include ../../Makefile.include

CC_FLAGS := -I../testing -I../cc $(CC_FLAGS)

ifeq ($(DEBUG),yes)
	CC_FLAGS := $(CC_FLAGS) -g
else
	CC_FLAGS := $(CC_FLAGS) -O2
endif

ifeq ($(OS),Windows_NT)
	CC_FLAGS := -mno-cygwin -mthreads $(CC_FLAGS)
endif


objects := $(patsubst %.cpp,%.o,$(wildcard *.cpp)) 
headers := $(wildcard *.h)

all: $(objects)

sylvester.a: $(objects)
	ar cr sylvester.a $(objects)
	ranlib sylvester.a

clear:
	rm -f *.o
	rm -f sylvester.a

%.o : %.cpp $(headers)
	$(CC) $(CC_FLAGS) $(EXTERN_DEFS) -c $*.cpp

