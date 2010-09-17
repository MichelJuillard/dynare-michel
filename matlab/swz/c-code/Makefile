# -*- Makefile -*- 
# Generic makefile for IMSL, MKL, blas/lapack, and FORTRAN libraries
# Generic makefile Tao's and Daniel's C routines 
# Generic makefile for various machines
# All output is to OUT_DIR

# MACHINE
STANDALONE_HB = STANDALONE_HB
#MACOSX = MACOSX

# LIBRARIES
USE_TAO = USE_TAO
#USE_MKL = USE_MKL
#USE_IMSL = USE_IMSL
USE_LAPACK = USE_LAPACK
USE_FORTRAN = USE_FORTRAN
USE_MATH = USE_MATH
USE_GSL = USE_GSL
#USE_GSL_CBLAS = USE_GSL_CBLAS


WORK_DIR = .
OUT_DIR = ../bin

CFLAGS = -g -DINTELCMATHLIBRARY -ansi

# MACHINE DEPENDENT DIRECTORIES 

# LAPTOP
ifdef USE_LAPTOP
WORK_DIR =/home/kalagman/Documents/InternshipJob/Dynare/
OUT_DIR = /home/kalagman/Documents/InternshipJob/Dynare/work/

CC = /usr/bin/gcc-4.3.2
CXX = g++
FC = gfortran
CFLAGS := $(CFLAGS) -static

TAO_DIR = $(WORK_DIR)/utilities/TZCcode
IMSL_BASE_DIR = /opt/vni/CTT6.0

GSL_INCLUDE_DIR = /usr/include/gsl
GSL_LIBS_DIR = /usr/lib

MKL_BASE_DIR = /usr
MKL_LIBS_DIR = $(MKL_BASE_DIR)/lib
MKL_LIBS = -llapack  -lblas

NPSOL_DIR = /home/f1dfw14/data/npsol-5.0

FORTRAN_LIBS = /usr/lib/gcc/i586-manbo-linux-gnu/4.3.2/libgfortran.a
endif

ifdef USE_MICHEL_LAPTOP

CC = gcc
CXX = g++
#CC = /opt/intel/Compiler/11.0/074/bin/intel64/icc -no-multibyte-chars
FC = gfortran
CFLAGS := $(CFLAGS) -static

TAO_DIR = $(WORK_DIR)/utilities/TZCcode

GSL_INCLUDE_DIR = /usr/include/gsl
GSL_LIBS_DIR = /usr/lib

INTEL_LIBS_DIR = /opt/intel/Compiler/11.0/074/lib/intel64
INTEL_LIBS = -lguide
MKL_BASE_DIR = /opt/intel/Compiler/11.0/074/mkl
MKL_LIBS_DIR = $(MKL_BASE_DIR)/lib/em64t
MKL_LIBS = -lmkl_lapack -lmkl_em64t

endif


ifdef STANDALONE_HB
CC = gcc
CXX = g++
FC = gfortran

GSL_INCLUDE_DIR = /usr/local64/include/gsl
GSL_LIBS_DIR = /usr/local64/lib

TAO_DIR = $(WORK_DIR)/utilities/TZCcode
endif

#ifdef MACOSX
#MAC_LIBS = #-framework Accelerate #-lmx -lmex -lmat -lstdc++ -L/Applications/MATLAB_R2009b.app/bin/maci64/ -lmwblas -lmwlapack 
#CFLAGS := $(CFLAGS) -DMACOSX
#endif


#################################################################################
VPATH = 
INCLUDE_DIR = -I/Applications/MATLAB_R2010b.app/extern/include
LIBS_DIR = 
LIBS = $(MAC_LIBS)
OBJS =
#################################################################################

### TAO'S FILES #################################################################
# TAO_DIR MUST be defined
ifdef USE_TAO
  VPATH := $(VPATH) $(TAO_DIR)
  INCLUDE_DIR := $(INCLUDE_DIR) -I$(TAO_DIR)
endif
#################################################################################

### GLS #########################################################################
# GSL_LIBS_DIR and GSL_INCLUDE_DIR muste be defined
ifdef USE_GSL
  INCLUDE_DIR := $(INCLUDE_DIR) -I$(GSL_INCLUDE_DIR)
  LIBS_DIR := $(LIBS_DIR) -L$(GSL_LIBS_DIR)
  LIBS := $(LIBS) -lgsl
endif
#################################################################################

### GLS CBLAS ###################################################################
# GSL_LIBS_DIR and GSL_INCLUDE_DIR muste be defined
ifdef USE_GSL_CBLAS
  ifndef USE_GSL
    INCLUDE_DIR := $(INCLUDE_DIR) -I$(GSL_INCLUDE_DIR)
    LIBS_DIR := $(LIBS_DIR) -L$(GSL_LIBS_DIR)
  endif
  LIBS := $(LIBS) -lgslcblas
endif
#################################################################################

### IMSL ########################################################################
# IMSL_BASE_DIR must be defined
ifdef USE_IMSL
  USE_PTHREAD = USE_PTHREAD 
  IMSL_LIBS_DIR =$(IMSL_BASE_DIR)/lib/lib.linux_icc
  INCLUDE_DIR := $(INCLUDE_DIR) -I$(IMSL_BASE_DIR)/include
  LIBS_DIR := $(LIBS_DIR) -L$(IMSL_LIBS_DIR)
  LIBS := $(LIBS) $(IMSL_LIBS_DIR)/libimslcmath.a $(IMSL_LIBS_DIR)/libimslcstat.a 
endif
#################################################################################

### NPSOL #######################################################################
ifdef USE_NPSOL
  NPSOL_LIBS = -lnpsol_f77 -llssol_f77
  LIBS_DIR := $(LIBS_DIR) -L$(NPSOL_DIR)
  LIBS := $(LIBS) $(NPSOL_LIBS)
endif
#################################################################################

### MKL #########################################################################
# MKL_BASE_DIR and MKL_LIBS must be defined
ifdef USE_MKL
  USE_PTHREAD = USE_PTHREAD 
  INCLUDE_DIR := $(INCLUDE_DIR) -I$(MKL_BASE_DIR)/include
  LIBS_DIR := $(LIBS_DIR) -Wl,-rpath=$(MKL_LIBS_DIR) -L$(MKL_LIBS_DIR) -L $(INTEL_LIBS_DIR)
  LIBS := $(LIBS) $(MKL_LIBS) $(INTEL_LIBS)
endif
#################################################################################

### BLAS / LAPACK LIBRARIES #####################################################
ifdef USE_LAPACK
  ifndef USE_MKL
    USE_FORTRAN = USE_FORTRAN
    ifdef LAPACK_LIBS_DIR
      LIBS_DIR := $(LIBS_DIR) -L$(LAPACK_LIBS_DIR)
    endif
    LIBS := $(LIBS) -llapack -lblas
  endif
endif
#################################################################################

### FORTRAN #####################################################################
# FORTRAN_LIBS should be defined
ifdef USE_FORTRAN
  ifdef FORTRAN_LIBS_DIR
    LIBS_DIR := $(LIBS_DIR) -L$(FORTRAN_LIBS_DIR)
  endif
  LIBS := $(LIBS) $(FORTRAN_LIBS)
endif
#################################################################################

### PTHREAD #####################################################################
ifdef USE_PTHREAD
  LIBS := $(LIBS) -lpthread
endif
#################################################################################

### MATH LIBRARY ################################################################
ifdef USE_MATH
  LIBS := $(LIBS) -lm
endif
#################################################################################

### DW FILES ####################################################################
MATRIX_DIR = $(WORK_DIR)/utilities/DWCcode/matrix
ERROR_DIR = $(WORK_DIR)/utilities/DWCcode/error
ARRAY_DIR = $(WORK_DIR)/utilities/DWCcode/arrays
ASCII_DIR = $(WORK_DIR)/utilities/DWCcode/ascii 
HIST_DIR = $(WORK_DIR)/utilities/DWCcode/histogram
STAT_DIR = $(WORK_DIR)/utilities/DWCcode/stat
SPHERICAL_DIR = $(WORK_DIR)/utilities/DWCcode/spherical
SORT_DIR = $(WORK_DIR)/utilities/DWCcode/sort
SWITCH_DIR = $(WORK_DIR)/sbvar/switching
VAR_DIR = $(WORK_DIR)/sbvar/var
#################################################################################

# DW  FILES
INCLUDE_DIR := -I$(MATRIX_DIR) -I$(ERROR_DIR) -I$(ARRAY_DIR) -I$(ASCII_DIR) -I$(HIST_DIR) -I$(STAT_DIR) -I$(SPHERICAL_DIR) -I$(SORT_DIR) -I$(SWITCH_DIR) -I$(VAR_DIR) $(INCLUDE_DIR) 
VPATH := $(VPATH) $(MATRIX_DIR) $(ERROR_DIR) $(ARRAY_DIR) $(ASCII_DIR) $(HIST_DIR) $(STAT_DIR) $(SPHERICAL_DIR) $(SORT_DIR) $(SWITCH_DIR) $(VAR_DIR)
OBJS := $(OBJS) bmatrix.o swzmatrix.o dw_error.o dw_rand.o dw_matrix_rand.o dw_array.o dw_matrix_array.o dw_matrix_sort.o dw_ascii.o dw_parse_cmd.o dw_histogram.o

# MEX
INCLUDE_DIR := -I$(WORK_DIR)/mex $(INCLUDE_DIR)
VPATH := $(VPATH) $(WORK_DIR)/mex
OBJS := $(OBJS) modify_for_mex.o

# TAO FILES
OBJS := $(OBJS)

# PROJECT FILE
INCLUDE_DIR := $(INCLUDE_DIR)
VPATH := $(VPATH)
OBJS_DRAWS := $(OBJS) PrintDraws.o switch.o switchio.o VARbase.o VARio.o command_line_VAR.o
OBJS_TAO := tzmatlab.o mathlib.o cstz_dw.o
OBJS_ESTIM := $(OBJS) $(OBJS_TAO) estimate.o VARbase.o VARio.o command_line_VAR.o switch.o switchio.o switch_opt.o csminwel.o
OBJS_INIT := $(OBJS) create_init_file.o switch.o switchio.o VARbase.o VARio.o VARio_matlab.o 
OBJS_MHM_1 := $(OBJS)  mhm_VAR_main_1.o mhm_VAR.o VARbase.o VARio.o command_line_VAR.o switch.o switchio.o 
OBJS_MHM_2 := $(OBJS) mhm_VAR_main_2.o spherical.o VARbase.o VARio.o switch.o switchio.o mhm_VAR.o
OBJS_PROBA := $(OBJS) probabilities.o switch.o switchio.o VARbase.o VARio.o command_line_VAR.o
OBJS_FORECAST := $(OBJS) forecast.o switch.o switchio.o VARbase.o VARio.o command_line_VAR.o


#################################################################################
# OUTPUT

all: $(OUT_DIR)/sbvar_draws $(OUT_DIR)/sbvar_estimation $(OUT_DIR)/sbvar_init_file $(OUT_DIR)/sbvar_mhm_1 $(OUT_DIR)/sbvar_mhm_2 $(OUT_DIR)/sbvar_probabilities $(OUT_DIR)/sbvar_forecast


$(OUT_DIR)/sbvar_draws: $(OBJS_DRAWS)
	$(CXX) $(CFLAGS) $^ $(LIBS_DIR) $(LIBS) -o $(OUT_DIR)/sbvar_draws

$(OUT_DIR)/sbvar_estimation: $(OBJS_ESTIM)
	$(CXX) $(CFLAGS) $^ $(LIBS_DIR) $(LIBS) -o $(OUT_DIR)/sbvar_estimation

$(OUT_DIR)/sbvar_init_file: $(OBJS_INIT)
	$(CXX) $(CFLAGS) $^ $(LIBS_DIR) $(LIBS) -o $(OUT_DIR)/sbvar_init_file

$(OUT_DIR)/sbvar_mhm_1: $(OBJS_MHM_1)
	$(CXX) $(CFLAGS) $^ $(LIBS_DIR) $(LIBS) -o $(OUT_DIR)/sbvar_mhm_1

$(OUT_DIR)/sbvar_mhm_2: $(OBJS_MHM_2)
	$(CXX) $(CFLAGS) $^ $(LIBS_DIR) $(LIBS) -o $(OUT_DIR)/sbvar_mhm_2

$(OUT_DIR)/sbvar_probabilities: $(OBJS_PROBA)
	$(CXX) $(CFLAGS) $^ $(LIBS_DIR) $(LIBS) -o $(OUT_DIR)/sbvar_probabilities

$(OUT_DIR)/sbvar_forecast: $(OBJS_FORECAST)
	$(CXX) $(CFLAGS) $^ $(LIBS_DIR) $(LIBS) -o $(OUT_DIR)/sbvar_forecast

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDE_DIR) -fexceptions -c $< -o $@

%.o : %.cc
	$(CXX) $(CFLAGS) $(INCLUDE_DIR) -fexceptions -c $< -o $@

#
#################################################################################
clean:
	rm -f $(OUT_DIR)/*.o *.o ; rm -f ../bin/*
#
#################################################################################
