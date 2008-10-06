# $Header: /var/lib/cvs/dynare_cpp/sylv/testing/Makefile,v 1.2 2004/09/28 16:16:43 kamenik Exp $

# Tag $Name:  $

LD_LIBS := -llapack -lcblas -lf77blas -latlas -lg2c
CC_FLAGS := -Wall -I../cc $(CC_INCLUDE_PATH)
CC_DEFS := #-DUSE_MEMORY_POOL
ifeq ($(DEBUG),yes)
	CC_FLAGS := $(CC_FLAGS) -g
	CC_DEFS := $(CC_DEFS) -DDEBUG
else
	CC_FLAGS := $(CC_FLAGS) -O2
endif
 
objects := $(patsubst %.cpp,%.o,$(wildcard ../cc/*.cpp)) 
headers := $(wildcard ../cc/*.h)
 
clear:
	make -C ../cc clear
	rm -f tests.exe
	rm -f tests
	rm -f *.o

../cc/%.o : ../cc/%.cpp $(headers)
	make EXTERN_DEFS="$(CC_DEFS)" -C ../cc $*.o

%.o : %.cpp $(headers)
	$(CC) $(CC_FLAGS) $(CC_DEFS) -c $*.cpp 

# untar testing data
tdata.done : tdata.tgz
	tar -xzf tdata.tgz
	touch tdata.done

test: $(objects) tests.o MMMatrix.o tdata.done
	$(CC) $(CC_FLAGS) -o test $(objects) tests.o MMMatrix.o $(LD_LIBS)
