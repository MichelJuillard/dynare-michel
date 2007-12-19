# $Header: /var/lib/cvs/dynare_cpp/sylv/matlab/Makefile,v 1.1.1.1 2004/06/04 13:01:13 kamenik Exp $

# Tag $Name:  $

# set directory to sylvester code
sylv=../cc

LD_LIBS := -llapack -lcblas -lblas -latlas
CC_FLAGS := -Wall -I../cc -fPIC
CC_DEFS := -I/usr/local/matlab73/extern/include -DMATLAB #-DUSE_MEMORY_POOL
MEX := mex.bat

objects := $(patsubst %.cpp,%.o,$(wildcard ../cc/*.cpp)) 
headers := $(wildcard ../cc/*.h)

# set mex file suffix
mex_suffix=dll
ifeq ($(OSTYPE),linux-gnu)
	mex_suffix=mexglx
endif

aa_gensylv.$(mex_suffix): gensylv.cpp ../cc/sylvester.a
	$(MEX) -I../cc -DMATLAB gensylv.cpp ../cc/sylvester.a c:/matlab6p5/extern/lib/win32/microsoft/msvc60/libmwlapack.lib
	mv gensylv.$(mex_suffix) aa_gensylv.$(mex_suffix)

../cc/sylvester.a : $(objects)
	make -C ../cc sylvester.a

../cc/%.o: ../cc/%.cpp $(headers)
	make EXTERN_DEFS="$(CC_DEFS)" -C ../cc $*.o

clear:
	make -C ../cc clear
	rm -f *.dll
	rm -f *.mexglx
