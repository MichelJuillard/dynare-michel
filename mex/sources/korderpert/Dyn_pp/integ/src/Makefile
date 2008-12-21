CC_FLAGS = -Wall -I../.. -I../../sylv/cc -I../../integ/cc -I../../tl/cc \
        $(CC_INCLUDE_PATH)
LDFLAGS = -llapack -lcblas -lf77blas -latlas -lg2c -lstdc++

ifeq ($(CC),)
	CC := gcc
endif

ifneq ($(LD_LIBRARY_PATH),)	# use LD_LIBRARY_PATH from environment
	LDFLAGS := -Wl,--library-path $(LD_LIBRARY_PATH) $(LDFLAGS)
endif

ifeq ($(DEBUG),yes)
	CC_FLAGS := $(CC_FLAGS) -g -DTL_DEBUG=2
else
	CC_FLAGS := $(CC_FLAGS) -O3 -DPOSIX_THREADS
endif

ifeq ($(OS),Windows_NT)
	CC_FLAGS := -mno-cygwin -mthreads $(CC_FLAGS)
	LDFLAGS := -mno-cygwin -mthreads $(LDFLAGS) -lpthreadGC2
	ARCH := w32
else
	LDFLAGS := $(LDFLAGS) -lpthread
	ARCH := linux
endif


sylvcppsource := $(wildcard ../../sylv/cc/*.cpp)
sylvhsource := $(wildcard ../../sylv/cc/*.h)
sylvobjects := $(patsubst %.cpp, %.o, $(sylvcppsource))


tlcwebsource := $(wildcard ../../tl/cc/*.cweb)
tlcppsource := $(patsubst %.cweb,%.cpp,$(tlcwebsource))
tlhwebsource := $(wildcard ../../tl/cc/*.hweb)
tlhsource := $(patsubst %.hweb,%.h,$(tlhwebsource))
tlobjects := $(patsubst %.cweb,%.o,$(tlcwebsource))

integcwebsource := $(wildcard ../../integ/cc/*.cweb)
integcppsource := $(patsubst %.cweb,%.cpp,$(integcwebsource))
integhwebsource := $(wildcard ../../integ/cc/*.hweb)
integhsource := $(patsubst %.hweb,%.h,$(integhwebsource))
integobjects := $(patsubst %.cweb,%.o,$(integcwebsource))

parserhsource := $(wildcard ../../parser/cc/*.h)
parsercppsource := $(wildcard ../../parser/cc/*.cpp)

utilshsource := $(wildcard ../../utils/cc/*.h)
utilscppsource := $(wildcard ../../utils/cc/*.cpp)
utilsobjects := $(patsubst %.cpp,%.o,$(utilscppsource))

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

../../utils/cc/%.o: ../../utils/cc/%.cpp $(utilshsource)
	make -C ../../utils/cc $*.o

../../parser/cc/%.o: ../../parser/cc/%.cpp $(parserhsource)
	make -C ../../parser/cc $*.o

../../parser/cc/matrix_tab.o:
	make -C ../../parser/cc matrix_tab.o

../../parser/cc/matrix_ll.o:
	make -C ../../parser/cc matrix_ll.o

quadrature-points: quadrature-points.cpp $(sylvhsource) $(sylvobjects) $(integhsource) $(integobjects) $(parserhsource) $(utilshsource) $(tlhsource) $(tlobjects) $(utilsobjects)
	$(CC) $(CC_FLAGS) quadrature-points.cpp -o quadrature-points ../../integ/cc/quadrature.o ../../integ/cc/smolyak.o ../../integ/cc/product.o ../../integ/cc/vector_function.o ../../tl/cc/sthread.o ../../tl/cc/symmetry.o ../../tl/cc/equivalence.o ../../tl/cc/int_sequence.o ../../tl/cc/tl_static.o ../../tl/cc/permutation.o ../../parser/cc/matrix_parser.o ../../parser/cc/matrix_tab.o ../../parser/cc/matrix_ll.o ../../parser/cc/parser_exception.o ../../sylv/cc/GeneralMatrix.o ../../sylv/cc/Vector.o ../../sylv/cc/SymSchurDecomp.o ../../sylv/cc/SylvException.o ../../utils/cc/memory_file.o $(LDFLAGS)
 

clear:
	make -C ../../tl/cc/ clear
	make -C ../../integ/cc clear
	make -C ../../parser/cc clear
	make -C ../../utils/cc clear
	make -C ../../sylv/cc clear
	rm -rf quadrature-points
