# $Header: /var/lib/cvs/dynare_cpp/sylv/cc/Makefile,v 1.4 2005/01/18 21:28:26 kamenik Exp $

# Tag $Name:  $
FF=f95

CC_FLAGS := -Wall -Winline -fpic
#CC_FLAGS := -Wall -Winline -I../testing -I../cc -DMATLAB  \
#        $(CC_INCLUDE_PATH) -I$(MATLAB_PATH)/extern/include #-pg

LDFLAGS = -Wl,-L$(MATLAB_PATH)/extern/lib/win32/microsoft/ \
	-Wl,-llibmex -Wl,-llibmx  -Wl,-llibmwlapack -Wl,-llibdflapack \
	-lg2c -lmingw32  -lstdc++ 

LD_LIBS=$(LDFLAGS)

ifeq ($(DEBUG),yes)
#	CC_FLAGS := $(CC_FLAGS) -g -DTL_DEBUG=2
    CC_FLAGS := $(CC_FLAGS) -g -DPOSIX_THREADS
else
	CC_FLAGS := $(CC_FLAGS) -O3
endif

ifeq ($(OS),Windows_NT)
	CC_FLAGS := -mno-cygwin -mthreads $(CC_FLAGS)
endif


objects := $(patsubst %.f90,%.o,$(wildcard *.f90)) 

all: $(objects)

clear:
	rm -f *.o

%.o : %.f90 
	$(FF) $(CC_FLAGS) $(EXTERN_DEFS) -c $*.f90

