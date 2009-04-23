all:
	make -C preprocessor
	make -C doc

clean:
	make -C preprocessor clean
	make -C doc clean
	rm -f matlab/dynare_m matlab/dynare_m.exe
	rm -f mex/2007a/* mex/2007b/* mex/octave/*.mex
	rm -f windows/*.exe
	rm -f config.log config.status *~
	rm -rf autom4te.cache
	rm -f preprocessor/Makefile

DYNAREBASE=$(shell basename $(shell pwd))

# Rule for creating the source tarball
# The parent directory must be called dynare-4.x.y
# WARNING: this rule will make your SVN working copy unusable!
srctarball: clean
	find -name .svn | xargs rm -rf
	find -type f -name '*~' | xargs rm -f
	cd ..; tar cvzf $(DYNAREBASE).tgz --exclude=tests $(DYNAREBASE)
