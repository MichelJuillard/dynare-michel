check_PROGRAMS = test-dr testModelSolution

test_dr_SOURCES = ../libmat/Matrix.cc ../libmat/Vector.cc ../libmat/QRDecomposition.cc ../libmat/GeneralizedSchurDecomposition.cc ../libmat/LUSolver.cc ../DecisionRules.cc test-dr.cc
test_dr_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS)
test_dr_CPPFLAGS = -I.. -I../libmat -I../../

testModelSolution_SOURCES = ../libmat/Matrix.cc ../libmat/Vector.cc ../libmat/QRDecomposition.cc ../libmat/GeneralizedSchurDecomposition.cc ../libmat/LUSolver.cc ../utils/dynamic_dll.cc ../DecisionRules.cc ../ModelSolution.cc testModelSolution.cc
testModelSolution_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) $(FLIBS) $(LIBADD_DLOPEN)
testModelSolution_CPPFLAGS = -I.. -I../libmat -I../../ -I../utils -DMEXEXT=\"$(MEXEXT)\"

check-local:
	./test-dr
