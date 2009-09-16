SUBDIRS =

if HAVE_MKOCTFILE
SUBDIRS += mjdgges/octave kronecker/octave gensylv/octave bytecode/octave
endif

if HAVE_MEXOPTS
SUBDIRS += mjdgges/matlab kronecker/matlab gensylv/matlab bytecode/matlab
endif

EXTRA_DIST = \
	../sources/matlab_versions_compatibility.h \
	../sources/build_matlab.m \
	../sources/build_matlab_multithread.m \
	../sources/build_octave.m
