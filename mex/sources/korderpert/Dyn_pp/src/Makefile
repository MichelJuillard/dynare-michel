CC_FLAGS = -Wall -I../sylv/cc -I../tl/cc -I../kord -I../integ/cc -I.. \
	$(CC_INCLUDE_PATH)

LDFLAGS = -llapack -lcblas -lf77blas -latlas -lg2c -lstdc++ -lcygwin
DYNVERSION = 1.3.7
tmpdir = /tmp/aaabbb123321
pthreadGC2 = `which pthreadGC2.dll`
mingwm10 = `which mingwm10.dll`

ifeq ($(CC),)
	CC := gcc
endif

ifneq ($(LD_LIBRARY_PATH),)	# use LD_LIBRARY_PATH from environment
#	LDFLAGS := -Wl,--library-path $(LD_LIBRARY_PATH) $(LDFLAGS)
	LDFLAGS := -Wl,--library-path $(LD_LIBRARY_LIST) $(LDFLAGS)
#	LDFLAGS := $(LD_LIBRARY_LIST) $(LDFLAGS)
endif

ifeq ($(DEBUG),yes)
	CC_FLAGS := $(CC_FLAGS) -g -DTL_DEBUG=2
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

sylvcppsource := $(wildcard ../sylv/cc/*.cpp)
sylvhsource := $(wildcard ../sylv/cc/*.h)
sylvobjects := $(patsubst %.cpp, %.o, $(sylvcppsource))

tlcwebsource := $(wildcard ../tl/cc/*.cweb)
tlcppsource := $(patsubst %.cweb,%.cpp,$(tlcwebsource))
tlhwebsource := $(wildcard ../tl/cc/*.hweb)
tlhsource := $(patsubst %.hweb,%.h,$(tlhwebsource))
tlobjects := $(patsubst %.cweb,%.o,$(tlcwebsource))

kordcwebsource := $(wildcard ../kord/*.cweb)
kordcppsource := $(patsubst %.cweb,%.cpp,$(kordcwebsource))
kordhwebsource := $(wildcard ../kord/*.hweb)
kordhsource := $(patsubst %.hweb,%.h,$(kordhwebsource))
kordobjects := $(patsubst %.cweb,%.o,$(kordcwebsource))

integcwebsource := $(wildcard ../integ/cc/*.cweb)
integcppsource := $(patsubst %.cweb,%.cpp,$(integcwebsource))
integhwebsource := $(wildcard ../integ/cc/*.hweb)
integhsource := $(patsubst %.hweb,%.h,$(integhwebsource))
integobjects := $(patsubst %.cweb,%.o,$(integcwebsource))

parserhsource := $(wildcard ../parser/cc/*.h)
parsercppsource := $(wildcard ../parser/cc/*.cpp)

utilshsource := $(wildcard ../utils/cc/*.h)
utilscppsource := $(wildcard ../utils/cc/*.cpp)

cppsource := $(wildcard *.cpp) $(patsubst %.y,%_ll.cc,$(wildcard *.y)) $(patsubst %.lex,%_tab.cc,$(wildcard *.lex))
hsource := $(wildcard *.h) $(patsubst %.lex,%_tab.hh,$(wildcard *.lex))
objects := $(patsubst %.cpp,%.o,$(wildcard *.cpp)) $(patsubst %.y,%_ll.o,$(wildcard *.y)) $(patsubst %.lex,%_tab.o,$(wildcard *.lex))

all: dynare++

../tl/cc/dummy.ch:
	make -C ../tl/cc dummy.ch

../tl/cc/%.cpp: ../tl/cc/%.cweb ../tl/cc/dummy.ch
	make -C ../tl/cc $*.cpp

../tl/cc/%.h: ../tl/cc/%.hweb ../tl/cc/dummy.ch
	make -C ../tl/cc $*.h

../tl/cc/%.o: ../tl/cc/%.cpp $(tlhsource)
	make -C ../tl/cc $*.o

../integ/cc/dummy.ch:
	make -C ../integ/cc dummy.ch

../integ/cc/%.cpp: ../integ/cc/%.cweb ../integ/cc/dummy.ch
	make -C ../integ/cc $*.cpp

../integ/cc/%.h: ../integ/cc/%.hweb ../integ/cc/dummy.ch
	make -C ../integ/cc $*.h

../integ/cc/%.o: ../integ/cc/%.cpp $(integhsource) $(tlhsource)
	make -C ../integ/cc $*.o


../sylv/cc/%.o: ../sylv/cc/%.cpp $(sylvhsource)
	make -C ../sylv/cc $*.o

../kord/dummy.ch:
	make -C ../kord dummy.ch

../kord/%.cpp: ../kord/%.cweb ../kord/dummy.ch
	make -C ../kord $*.cpp

../kord/%.h: ../kord/%.hweb ../kord/dummy.ch
	make -C ../kord $*.h

../kord/%.o: ../kord/%.cpp $(tlhsource) $(kordhsource) $(integhsource)
	make -C ../kord $*.o

# for some reason, the pattern rules didn't work here, so the rules
# are expanded:
dynglob_tab.cc: dynglob.y dynare_atoms.h dynare_model.h
	bison -d -t --verbose -odynglob_tab.cc dynglob.y

dynglob_tab.hh: dynglob.y dynare_atoms.h dynare_model.h
	bison -d -t --verbose -odynglob_tab.cc dynglob.y

dynglob_ll.cc: dynglob.lex dynglob_tab.hh
	flex -i -odynglob_ll.cc dynglob.lex

dynglob_ll.o: dynglob_ll.cc $(hsource)
	$(CC) -I.. -I../sylv/cc -I../tl/cc -O3 -c dynglob_ll.cc

dynglob_tab.o: dynglob_tab.cc $(hsource)
	$(CC) -I.. -I../sylv/cc -I../tl/cc -O3 -c dynglob_tab.cc

%.o: %.cpp $(hsource) $(kordhsource) $(tlhsource) $(sylvhsource) $(parserhsource)
	$(CC) $(CC_FLAGS) -DDYNVERSION=\"$(DYNVERSION)\" -o $*.o -c $*.cpp

../parser/cc/parser.a: $(parserhsource) $(parsercppsource)
	make -C ../parser/cc parser.a

../utils/cc/utils.a: $(utilshsource) $(utilscppsource)
	make -C ../utils/cc utils.a

dynare++: $(tlhwebsource) $(tlcwebsoure) $(tlhsource) $(tlcppsource) \
         $(integhwebsource) $(integcwebsoure) $(integhsource) $(integcppsource) \
         $(kordhwebsource) $(kordcwebsoure) $(kordhsource) $(kordcppsource) \
         $(sylvhsource) $(sylvcppsource) \
         $(kordobjects) $(tlobjects) $(integobjects) $(sylvobjects) $(objects) \
         ../parser/cc/parser.a ../utils/cc/utils.a
	$(CC) -g -o dynare++ $(objects) $(kordobjects) $(integobjects) $(tlobjects) ../parser/cc/parser.a ../utils/cc/utils.a $(sylvobjects) $(LDFLAGS)

../extern/matlab/dynare_simul_.$(MEX_SUFFIX):
	make -C ../extern/matlab all

srcball:
	rm -rf $(tmpdir)
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/tests
	cp ../tests/*.mod ../tests/*.dyn $(tmpdir)/dynare++-$(DYNVERSION)/tests
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/src
	cp *.h *.cpp *.y *.lex Makefile $(tmpdir)/dynare++-$(DYNVERSION)/src
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/tl/cc
	cp ../tl/cc/*web ../tl/cc/Makefile $(tmpdir)/dynare++-$(DYNVERSION)/tl/cc
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/tl/testing
	cp ../tl/testing/*.cpp ../tl/testing/*.h ../tl/testing/Makefile $(tmpdir)/dynare++-$(DYNVERSION)/tl/testing
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/integ/cc
	cp ../integ/cc/*web ../integ/cc/Makefile ../integ/cc/*.dat $(tmpdir)/dynare++-$(DYNVERSION)/integ/cc
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/integ/testing
	cp ../integ/testing/*.cpp ../integ/testing/Makefile $(tmpdir)/dynare++-$(DYNVERSION)/integ/testing
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/sylv/cc
	cp ../sylv/cc/*.cpp ../sylv/cc/*.h ../sylv/cc/Makefile $(tmpdir)/dynare++-$(DYNVERSION)/sylv/cc
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/sylv/testing
	cp ../sylv/testing/*.cpp ../sylv/testing/*.h ../sylv/testing/Makefile $(tmpdir)/dynare++-$(DYNVERSION)/sylv/testing
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/kord
	cp ../kord/*web ../kord/Makefile ../kord/tests.cpp $(tmpdir)/dynare++-$(DYNVERSION)/kord
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/parser/cc
	cp ../parser/cc/*.cpp ../parser/cc/*.h ../parser/cc/*.lex ../parser/cc/*.y ../parser/cc/Makefile $(tmpdir)/dynare++-$(DYNVERSION)/parser/cc
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/utils/cc
	cp ../utils/cc/*.cpp ../utils/cc/*.h ../utils/cc/Makefile $(tmpdir)/dynare++-$(DYNVERSION)/utils/cc

	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/extern/matlab
	cp ../extern/matlab/*.cpp ../extern/matlab/*.m ../extern/matlab/Makefile ../extern/matlab/*.bat $(tmpdir)/dynare++-$(DYNVERSION)/extern/matlab
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/doc
	cp ../doc/*.tex $(tmpdir)/dynare++-$(DYNVERSION)/doc
	cp ../change_log.html $(tmpdir)/dynare++-$(DYNVERSION)/
	cd $(tmpdir); tar czf dynare++-$(DYNVERSION).src.tgz dynare++-$(DYNVERSION)/*
	mv $(tmpdir)/dynare++-$(DYNVERSION).src.tgz ..
	rm -rf $(tmpdir)

binball: dynare++ ../extern/matlab/dynare_simul_.$(MEX_SUFFIX)
	rm -rf $(tmpdir)
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/tests
	cp ../tests/*.mod ../tests/*.dyn $(tmpdir)/dynare++-$(DYNVERSION)/tests
	mkdir -p $(tmpdir)/dynare++-$(DYNVERSION)/doc
	cp ../doc/*.pdf $(tmpdir)/dynare++-$(DYNVERSION)/doc
	cp dynare++* $(tmpdir)/dynare++-$(DYNVERSION)
	cp ../extern/matlab/*.m $(tmpdir)/dynare++-$(DYNVERSION)
	cp ../extern/matlab/dynare_simul_.$(MEX_SUFFIX) $(tmpdir)/dynare++-$(DYNVERSION)
	cp ../change_log.html $(tmpdir)/dynare++-$(DYNVERSION)
ifeq ($(OS),Windows_NT)
	cp $(pthreadGC2) $(tmpdir)/dynare++-$(DYNVERSION)
	cp $(mingwm10) $(tmpdir)/dynare++-$(DYNVERSION)
	cd $(tmpdir); zip -r dynare++-$(DYNVERSION)-$(ARCH).ix86.zip dynare++-$(DYNVERSION)/*
	mv $(tmpdir)/dynare++-$(DYNVERSION)-$(ARCH).ix86.zip ..
else
	cd $(tmpdir); tar czf dynare++-$(DYNVERSION)-$(ARCH).ix86.tgz dynare++-$(DYNVERSION)/*
	mv $(tmpdir)/dynare++-$(DYNVERSION)-$(ARCH).ix86.tgz ..
endif
	rm -rf $(tmpdir)


clear:
	rm -f *.o
	rm -f dynare++*
	rm -f *_ll.cc *_tab.hh *_tab.cc *.output
	make -C ../tl/testing clear
	make -C ../tl/cc clear
	make -C ../integ/testing clear
	make -C ../integ/cc clear
	make -C ../sylv/testing clear
	make -C ../sylv/cc clear
	make -C ../kord clear
	make -C ../parser/cc clear
	make -C ../utils/cc clear
