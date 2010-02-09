check_PROGRAMS = test-dr

test_dr_SOURCES = ../libmat/Matrix.cc ../libmat/Vector.cc ../libmat/QRDecomposition.cc ../libmat/GeneralizedSchurDecomposition.cc ../libmat/LUSolver.cc ../DecisionRules.cc test-dr.cc
test_dr_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
test_dr_CPPFLAGS = -I.. -I../libmat -I../../
