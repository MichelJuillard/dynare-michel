SUBDIRS = preprocessor doc tests mex/sources

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
	contrib \
	NEWS \
	license.txt \
	windows/dynare.nsi \
	windows/mexopts-win32.bat \
	windows/mexopts-win64.bat \
	windows/README.txt \
	osx \
	examples \
	scripts

dist-hook:
	rm -rf `find $(distdir)/matlab $(distdir)/examples -name *~`
	rm -f $(distdir)/matlab/dynare_m$(EXEEXT) $(distdir)/matlab/dynare_version.m
	$(MKDIR_P) $(distdir)/mex/matlab $(distdir)/mex/octave
	rm -rf `find $(distdir)/contrib -name '.git*'`
