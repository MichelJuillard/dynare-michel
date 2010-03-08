check_PROGRAMS = test-qr test-gsd test-lu test-repmat

test_qr_SOURCES = ../Matrix.cc ../Vector.cc ../QRDecomposition.cc test-qr.cc
test_qr_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
test_qr_CPPFLAGS = -I.. -I../../../

test_gsd_SOURCES = ../Matrix.cc ../Vector.cc ../GeneralizedSchurDecomposition.cc test-gsd.cc
test_gsd_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
test_gsd_CPPFLAGS = -I.. -I../../../

test_lu_SOURCES = ../Matrix.cc ../Vector.cc ../LUSolver.cc test-lu.cc
test_lu_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
test_lu_CPPFLAGS = -I.. -I../../../

test_repmat_SOURCES = ../Matrix.cc ../Vector.cc test-repmat.cc
test_repmat_CPPFLAGS = -I..

check-local:
	./test-qr
	./test-gsd
	./test-lu
	./test-repmat
