# $Id: Makefile 2344 2009-02-09 20:36:08Z michel $
# Copyright 2004, Ondra Kamenik

include ../Makefile.include

CC_FLAGS := -I../tl/cc -I../sylv/cc -I../integ/cc $(CC_FLAGS)

xbsylvcppsource := $(wildcard ../sylv/cc/*.cpp)
sylvhsource := $(wildcard ../sylv/cc/*.h)
sylvobjects := $(patsubst %.cpp, %.o, $(sylvcppsource)) 
tlcwebsource := $(wildcard ../tl/cc/*.cweb)
tlcppsource := $(patsubst %.cweb,%.cpp,$(tlcwebsource))
tlhwebsource := $(wildcard ../tl/cc/*.hweb)
tlhsource := $(patsubst %.hweb,%.h,$(tlhwebsource))
tlobjects := $(patsubst %.cweb,%.o,$(tlcwebsource))
integhwebsource := $(wildcard ../integ/cc/*.hweb)
integhsource := $(patsubst %.hweb,%.h,$(integhwebsource))
cwebsource := $(wildcard *.cweb)
cppsource := $(patsubst %.cweb,%.cpp,$(cwebsource)) 
objects := $(patsubst %.cweb,%.o,$(cwebsource))
hwebsource := $(wildcard *.hweb)
hsource := $(patsubst %.hweb,%.h,$(hwebsource))

../integ/cc/dummy.ch:
	make -C ../integ/cc dummy.ch

../tl/cc/dummy.ch:
	make -C ../tl/cc dummy.ch

../tl/cc/%.cpp: ../tl/cc/%.cweb ../tl/cc/dummy.ch
	make -C ../tl/cc $*.cpp

../tl/cc/%.h: ../tl/cc/%.hweb ../tl/cc/dummy.ch
	make -C ../tl/cc $*.h

../integ/cc/%.h: ../integ/cc/%.hweb ../integ/cc/dummy.ch
	make -C ../integ/cc $*.h

../tl/cc/%.o: ../tl/cc/%.cpp $(tlhsource)
	make -C ../tl/cc $*.o

../sylv/cc/%.o: ../sylv/cc/%.cpp $(sylvhsource)
	make -C ../sylv/cc $*.o

dummy.ch:
	touch dummy.ch

%.cpp: %.cweb dummy.ch
	ctangle -bhp $*.cweb dummy.ch $*.cpp

%.h: %.hweb dummy.ch
	ctangle -bhp $*.hweb dummy.ch $*.h

%.o: %.cpp $(hwebsource) $(hsource) $(tlhwebsource) $(tlhsource) \
     $(integhwebsource) $(integhsource)
	$(CC) $(CC_FLAGS) -c $*.cpp

tests: $(hwebsource) $(cwebsoure) $(hsource) $(cppsource) \
       $(tlhwebsource) $(tlcwebsoure) $(tlhsource) $(tlcppsource) \
       $(sylvhsource) $(sylvcppsource) \
       tests.o $(objects) $(tlobjects) $(sylvobjects)
	$(CC) $(CC_FLAGS) $(objects) $(tlobjects) $(sylvobjects) tests.o -o tests $(LD_LIBS) 

kord.pdf: doc

doc: main.web $(hwebsource) $(cwebsource)
	cweave -bhp main.web
	pdftex main
	mv main.pdf kord.pdf 

clear:
	rm -f $(hsource)
	rm -f $(cppsource)
	rm -f *.o
	rm -f tests
	make -C ../tl/cc clear
	make -C ../sylv/cc clear
	rm -f main.{idx,dvi,pdf,scn,log,tex,toc}
