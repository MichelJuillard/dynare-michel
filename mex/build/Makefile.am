SUBDIRS =

if HAVE_MKOCTFILE
SUBDIRS += mjdgges/octave kronecker/octave gensylv/octave bytecode/octave
endif

if HAVE_MEXOPTS
SUBDIRS += mjdgges/matlab kronecker/matlab gensylv/matlab bytecode/matlab
endif

if HAVE_BLAS
if HAVE_LAPACK
if HAVE_PTHREAD
if HAVE_MEXOPTS
#SUBDIRS += korderpert/matlab korderpert/tests
endif
if HAVE_MKOCTFILE
#SUBDIRS += korderpert/octave
endif
endif
endif
endif


EXTRA_DIST = mex.def mexFunction-MacOSX.map
