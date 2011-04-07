# -*- Makefile -*- 
# Generic makefile for IMSL, MKL, blas/lapack, and FORTRAN libraries
# Generic makefile Tao's and Daniel's C routines 
# Generic makefile for various machines

# MACHINE
#USE_HOME=USE_HOME
#USE_LAPTOP=USE_LAPTOP 
USE_WORK_64=USE_WORK_64
#USE_WORK_32=USE_WORK_32
#USE_WULF=USE_WULF

# LIBRARIES SUPPORTED
#USE_TZ = USE_TZ
USE_DW = USE_DW
  USE_DW_MATRIX = USE_DW_MATRIX
  USE_DW_ERROR = USE_DW_ERROR
  USE_DW_ARRAY = USE_DW_ARRAY
  USE_DW_ASCII = USE_DW_ASCII
  USE_DW_STAT = USE_DW_STAT
  USE_DW_SORT = USE_DW_SORT
  USE_DW_ELLIPTICAL = USE_DW_ELLIPTICAL
  USE_DW_HISTO = USE_DW_HISTO
  USE_DW_MATH = USE_DW_MATH
  USE_DW_SWITCH = USE_DW_SWITCH
  USE_DW_STATE_SPACE = USE_DW_STATE_SPACE
#USE_MKL = USE_MKL
#USE_ATLAS = USE_ATLAS
#USE_IMSL = USE_IMSL
USE_LAPACK = USE_LAPACK
#USE_FORTRAN = USE_FORTRAN
USE_MATH = USE_MATH
USE_GSL = USE_GSL
#USE_GSL_CBLAS = USE_GSL_CBLAS
#USE_BOOST = USE_BOOST
#USE_MYSQL_CPP_CONN = USE_MYSQL_CPP_CONN
#USE_LOG4CXX = USE_LOG4CXX

BASE_DIR = .

USE_GCC = USE_GCC
CFLAGS = -g -DSTRUCTURED_COLUMN_MAJOR
#CFLAGS = -g -static

MAKEFILE_DIR = $(BASE_DIR)

include $(MAKEFILE_DIR)/make_machines_master
include $(MAKEFILE_DIR)/make_libraries_master

#################################################################################
# PROJECT PATHS
PROJECT_DIR = $(BASE_DIR)
PROJECT_SRC_DIR = $(BASE_DIR)/switch_dw/state_space/sbvar
EXE_DIR = $(BASE_DIR)
BUILD_DIR = $(BASE_DIR)

INCLUDE_DIR := $(INCLUDE_DIR) -I$(PROJECT_SRC_DIR)
VPATH := $(VPATH) $(PROJECT_SRC_DIR) $(BUILD_DIR)  

# INIT_EXEC #####################################################################
# DW FILES
INIT_OBJS = $(OBJS) dw_switch.o dw_switchio.o dw_dirichlet_restrictions.o dw_metropolis_theta.o

# PROJECT FILES
INIT_OBJS := $(INIT_OBJS) create_init_file.o VARbase.o VARio.o VARio_matlab.o

# OUTPUT
INIT_EFILE = $(EXE_DIR)/sbvar_init_file


# EXEC ######################################################################
# DW FILES
EXEC_OBJS = $(OBJS) dw_switch.o dw_switchio.o dw_dirichlet_restrictions.o dw_metropolis_theta.o dw_switch_opt.o dw_mdd_switch.o

# TZ FILES
# EXEC_OBJS := $(EXEC_OBJS)  tzmatlab.o csminwel.o

# PROJECT FILES
EXEC_OBJS := $(EXEC_OBJS) VARbase.o VARio.o dw_sbvar_command_line.o sbvar_estimate.o sbvar_simulate.o \
                          sbvar_probabilities.o sbvar_mdd.o dw_csminwel.o sbvar_impulse_responses.o sbvar_forecast.o

# OUTPUT
EXEC_EFILE = $(EXE_DIR)/sbvar_commandline


#################################################################################
all : INIT_EXEC EXEC

INIT_EXEC : $(INIT_EFILE)

EXEC : $(EXEC_EFILE)

$(INIT_EFILE): $(INIT_OBJS)
	$(CC) $(CFLAGS) $^ $(LIBS_DIR) $(LIBS) -o $(INIT_EFILE)

$(EXEC_EFILE): $(EXEC_OBJS)
	$(CC) $(CFLAGS) $^ $(LIBS_DIR) $(LIBS) -o $(EXEC_EFILE)

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDE_DIR) -c $< -o $(BUILD_DIR)/$@

%.o : %.cpp
	$(CC) $(CFLAGS) $(INCLUDE_DIR) -c $< -o $(BUILD_DIR)/$@

#################################################################################
clean:
	rm -f $(BUILD_DIR)/*.o sbvar_commandline sbvar_init_file *~
#
#################################################################################
test:
	echo "EXE_DIR: $(EXE_DIR)"
	echo "BUILD_DIR: $(BUILD_DIR)"
	echo "VPATH: $(VPATH)"
	echo "INCLUDE_DIR: $(INCLUDE_DIR)"
	echo "LIBS_DIR: $(LIBS_DIR)"
	echo "OBJS: $(OBJS)"
	echo "INIT_OBJS: $(INIT_OBJS)"
	echo "EXEC_OBJS: $(EXEC_OBJS)"
