CC = gcc
DEBUG = yes
MATLAB =1
#CC_INCLUDE_PATH=$CC_INCLUDE_PATH -I'c:/"Program Files"/MATLAB_SV71/extern/include'
CC_FLAGS = -Wall -I../../sylv/cc -I../../tl/cc -I../../kord -I../../integ/cc \
	 $(CC_INCLUDE_PATH) -Ic:/"Program Files"/MATLAB_SV71/extern/include

#LDFLAGS = -llapack -lcblas -lf77blas -latlas -lg2c -lstdc++

LDFLAGS = -Wl,-L"c:/Program Files"/MATLAB_SV71/extern/lib/win32/microsoft/ \
	-Wl,-llibmex -Wl,-llibmx  -Wl,-llibmwlapack -Wl,-llibdflapack \
	-lg2c -lmingw32  -lstdc++ 

ifeq ($(DEBUG),yes)
#	CC_FLAGS := $(CC_FLAGS) -g -DTL_DEBUG=2
    CC_FLAGS := $(CC_FLAGS) -g -DMATLAB -DPOSIX_THREADS
else
	CC_FLAGS := $(CC_FLAGS) -O3 -DPOSIX_THREADS
endif

ifeq ($(OS),Windows_NT)
	CC_FLAGS := -mno-cygwin -mthreads $(CC_FLAGS)
	LDFLAGS := -mno-cygwin -mthreads $(LDFLAGS)  -lpthreadGC2
	ARCH := w32
	MEX_SUFFIX = dll
else
	LDFLAGS := $(LDFLAGS) -lpthread
	ARCH := linux
	MEX_SUFFIX = mexglx
endif

sylvcppsource := $(wildcard ../../sylv/cc/*.cpp)
sylvhsource := $(wildcard ../../sylv/cc/*.h)
sylvobjects := $(patsubst %.cpp, %.o, $(sylvcppsource))

tlcppsource := $(wildcard ../../tl/cc/*.cpp)
tlhsource := $(wildcard ../../tl/cc/*.h)
tlobjects := $(patsubst %.cpp,%.o,$(tlcppsource))

#kordcwebsource := $(wildcard ../../kord/*.cweb)
kordcppsource := $(wildcard ../../kord/*.cpp) # $(patsubst %.cweb,%.cpp,$(kordcwebsource))
#kordhwebsource := $(wildcard ../../kord/*.hweb)
kordhsource := $(wildcard ../../kord/*.h)# $(patsubst %.hweb,%.h,$(kordhwebsource))
kordobjects := $(patsubst %.cpp,%.o,$(kordcppsource))

#integcwebsource := $(wildcard ../../integ/cc/*.cweb)
#integcppsource := $(patsubst %.cweb,%.cpp,$(integcwebsource))
#integhwebsource := $(wildcard ../../integ/cc/*.hweb)
#integhsource := $(patsubst %.hweb,%.h,$(integhwebsource))
#integobjects := $(patsubst %.cweb,%.o,$(integcwebsource))

cppsource := $(wildcard *.cpp)
mexobjects := $(patsubst %.cpp,%_.$(MEX_SUFFIX),$(cppsource))

all: $(mexobjects)

../../tl/cc/dummy.ch:
	make -C ../../tl/cc dummy.ch

../../tl/cc/%.cpp: ../../tl/cc/%.cweb ../../tl/cc/dummy.ch
	make -C ../../tl/cc $*.cpp

../../tl/cc/%.h: ../../tl/cc/%.hweb ../../tl/cc/dummy.ch
	make -C ../../tl/cc $*.h

../../tl/cc/%.o: ../../tl/cc/%.cpp $(tlhsource)
	make -C ../../tl/cc $*.o

../../integ/cc/dummy.ch:
	make -C ../../integ/cc dummy.ch

../../integ/cc/%.cpp: ../../integ/cc/%.cweb ../../integ/cc/dummy.ch
	make -C ../../integ/cc $*.cpp

../../integ/cc/%.h: ../../integ/cc/%.hweb ../../integ/cc/dummy.ch
	make -C ../../integ/cc $*.h

../../integ/cc/%.o: ../../integ/cc/%.cpp $(integhsource) $(tlhsource)
	make -C ../../integ/cc $*.o


../../sylv/cc/%.o: ../../sylv/cc/%.cpp $(sylvhsource)
	make -C ../../sylv/cc $*.o

../../kord/dummy.ch:
	make -C ../../kord dummy.ch

../../kord/%.cpp: ../../kord/%.cweb ../../kord/dummy.ch
	make -C ../../kord $*.cpp

../../kord/%.h: ../../kord/%.hweb ../../kord/dummy.ch
	make -C ../../kord $*.h

../../kord/%.o: ../../kord/%.cpp $(tlhsource) $(kordhsource) $(integhsource)
	make -C ../../kord $*.o

k_ord_lib.a: $(tlhsource) $(tlcppsource) \
          $(integhsource) $(integcppsource) \
          $(kordhsource) $(kordcppsource) \
         $(sylvhsource) $(sylvcppsource) \
         $(kordobjects) $(tlobjects) $(integobjects) $(sylvobjects)
	ar cr k_ord_lib.a  $(kordobjects) $(tlobjects) $(integobjects) $(sylvobjects)
	ranlib k_ord_lib.a


# to compile mex objects for Windows do:
# 1. install gnumex
# 2. create mexopts.bat via gnumex in this directory, specify MinGW compilation, and dll output
# 3. change -Ic:/MATLAB71/extern/include according your Matlab setup
# 4. pray it works
# OR: just use the mexopt.bat from the repository and check MATLAB root directory
%_.$(MEX_SUFFIX): %.cpp  $(tlhsource) $(tlcppsource) \
          $(integhsource) $(integcppsource) \
          $(kordhsource) $(kordcppsource) \
         $(sylvhsource) $(sylvcppsource) \
         k_ord_lib.a.a
ifeq ($(OS),Windows_NT)
	mex.bat -I../../sylv/cc -I../../tl/cc -I../../kord -I../../integ/cc COMPFLAGS\#"-c -DMATLAB_MEX_FILE -DPOSIX_THREADS" OPTIMFLAGS\#"-O3 -fexceptions -Ic:/MATLAB71/extern/include" GM_ADD_LIBS\#"$(LDFLAGS)" $*.cpp k_ord_lib.a.a
else
	mex -I../../sylv/cc/ -I../../tl/cc -I../../kord -I../../integ/cc $*.cpp CFLAGS='$$CFLAGS -fexceptions' k_ord_lib.a.a
endif
	mv $*.$(MEX_SUFFIX) $*_.$(MEX_SUFFIX)

clear:
	rm -f k_ord_lib.a.a
	rm -f *.mexglx
	rm -f *.dll
	make -C ../../tl/testing clear
	make -C ../../tl/cc clear
	make -C ../../integ/testing clear
	make -C ../../integ/cc clear
	make -C ../../sylv/testing clear
	make -C ../../sylv/cc clear
	make -C ../../kord clear
