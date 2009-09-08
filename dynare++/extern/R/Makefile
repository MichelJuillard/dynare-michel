RINTERNALS=/usr/share/R/include/

sylvcppsource := $(wildcard ../../sylv/cc/*.cpp)
sylvhsource := $(wildcard ../../sylv/cc/*.h)
sylvobjects := $(patsubst %.cpp, %.o, $(sylvcppsource))

tlcwebsource := $(wildcard ../../tl/cc/*.cweb)
tlcppsource := $(patsubst %.cweb,%.cpp,$(tlcwebsource))
tlhwebsource := $(wildcard ../../tl/cc/*.hweb)
tlhsource := $(patsubst %.hweb,%.h,$(tlhwebsource))
tlobjects := $(patsubst %.cweb,%.o,$(tlcwebsource))

kordcwebsource := $(wildcard ../../kord/*.cweb)
kordcppsource := $(patsubst %.cweb,%.cpp,$(kordcwebsource))
kordhwebsource := $(wildcard ../../kord/*.hweb)
kordhsource := $(patsubst %.hweb,%.h,$(kordhwebsource))
kordobjects := $(patsubst %.cweb,%.o,$(kordcwebsource))

integcwebsource := $(wildcard ../../integ/cc/*.cweb)
integcppsource := $(patsubst %.cweb,%.cpp,$(integcwebsource))
integhwebsource := $(wildcard ../../integ/cc/*.hweb)
integhsource := $(patsubst %.hweb,%.h,$(integhwebsource))
integobjects := $(patsubst %.cweb,%.o,$(integcwebsource))

parserhsource := $(wildcard ../../parser/cc/*.h)
parsercppsource := $(wildcard ../parser/cc/*.cpp)

utilshsource := $(wildcard ../../utils/cc/*.h)
utilscppsource := $(wildcard ../utils/cc/*.cpp)

srccpp := dynare3.cpp dynare_model.cpp planner_builder.cpp dynare_atoms.cpp dynare_params.cpp  nlsolve.cpp
objects := $(patsubst %.cpp,../../src/%.o,$(srccpp)) \
$(patsubst %.y,%_ll.o,$(wildcard ../../src/*.y)) \
$(patsubst %.lex,%_tab.o,$(wildcard ../../src/*.lex))

PKG_CPPFLAGS= -I../../tl/cc -I../../sylv/cc -I../../kord -I../../src -I../.. -I$(RINTERNALS)
PKG_LIBS= ${LAPACK_LIBS} ${BLAS_LIBS} $(objects) $(kordobjects) $(integobjects) $(tlobjects) ../../parser/cc/parser.a ../../utils/cc/utils.a $(sylvobjects) -lpthread -llapack -lcblas -lf77blas -latlas -lg2c -lstdc++

ifneq ($(LD_LIBRARY_PATH),)	# use LD_LIBRARY_PATH from environment
	PKG_LIBS := -Wl,--library-path $(LD_LIBRARY_PATH) $(PKG_LIBS)
endif

dynareR.so: dynareR.o
	g++ -shared  -o dynareR.so dynareR.o -L/usr/lib/R/lib -lR $(PKG_LIBS)

dynareR.o: dynareR.cpp
	g++ -I/usr/share/R/include -I/usr/share/R/include $(PKG_CPPFLAGS) \
	-fpic  -g -O2 -c dynareR.cpp -o dynareR.o -DDEBUG

test: test.cpp dynareR.cpp
	g++ -O0 -g -o test test.cpp -DDEBUG $(PKG_LIBS) $(PKG_CPPFLAGS)

test-debug:
	valgrind --leak-check=yes ./test
