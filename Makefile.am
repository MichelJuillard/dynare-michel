SUBDIRS = preprocessor doc tests

if HAVE_BLAS
if HAVE_LAPACK
SUBDIRS += dynare++
endif
endif

# MEX must be built after dynare++ (because of kordepert)
if ENABLE_MATLAB
SUBDIRS += mex/build/matlab
endif
if ENABLE_OCTAVE
SUBDIRS += mex/build/octave
endif

ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = \
	matlab \
	mex/sources \
	uncrustify.cfg \
	windows/dynare.nsi \
	windows/mexopts.bat \
	windows/README.txt

clean-local:
	rm -rf `find mex/sources -name *.o`

dist-hook: clean-local
	rm -rf `find $(distdir) -name .svn`
	rm -rf `find $(distdir)/mex/sources -name *.o`
	rm -f $(distdir)/matlab/dynare_m$(EXEEXT)
	$(MKDIR_P) $(distdir)/mex/matlab $(distdir)/mex/octave
