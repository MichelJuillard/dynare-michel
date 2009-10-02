SUBDIRS = preprocessor doc tests

if HAVE_BLAS
if HAVE_LAPACK
SUBDIRS += dynare++
endif
endif

# MEX must be built after dynare++ (because of kordepert)
SUBDIRS += mex/build/octave mex/build/matlab

ACLOCAL_AMFLAGS = -I m4

EXTRA_DIST = \
	matlab \
	mex/sources \
	uncrustify.cfg \
	windows/dynare.nsi \
	windows/fdl.txt \
	windows/gpl.txt \
	windows/README.txt

clean-local:
	rm -rf `find mex/sources -name *.o`

dist-hook: clean-local
	rm -rf `find $(distdir) -name .svn`
	rm -rf `find $(distdir)/mex/sources -name *.o`
	rm -f $(distdir)/matlab/dynare_m$(EXEEXT)
	$(MKDIR_P) $(distdir)/mex/2007a $(distdir)/mex/2007b $(distdir)/mex/2009a-64bit $(distdir)/mex/octave
