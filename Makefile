# Makefile for creating the source tarball
# The parent directory must be called dynare-4.x.y
# WARNING: this makefile will destroy all your .svn subdirectories!

DYNAREBASE=$(shell basename $(shell pwd))

srctarball:
	make -C preprocessor clean
	make -C doc clean
	rm -f matlab/dynare_m matlab/dynare_m.exe
	rm -f mex/2007a/* mex/2007b/* mex/octave/*.mex
	rm -f windows/*.exe
	find -name .svn | xargs rm -rf
	find -type f -name '*~' | xargs rm -f
	cd ..; tar cvzf $(DYNAREBASE).tgz \
		$(DYNAREBASE)/preprocessor \
	  $(DYNAREBASE)/matlab \
    $(DYNAREBASE)/doc \
    $(DYNAREBASE)/mex
