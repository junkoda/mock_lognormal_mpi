#
# mock_lognormal_mpi
#   

all: mock_lognormal_mpi

#
# Compiler configurations
#
CXX = mpic++

OPT := -DDOUBLEPRECISION

# Define OPENMP to enable MPI+OpenMP hybrid parallelization

# CXX is defined in ../Makefile

WOPT    ?= -Wall
CPPFLAGS  := -O3 $(WOPT) $(OPENMP) $(OPT) -fPIC
LIBS    := -lm 

# Define paths of FFTW3 & GSL libraries if necessary.

FFTW3_DIR ?= #e.g. /Users/jkoda/Research/opt/gcc/fftw3
GSL_DIR   ?= #e.g. /Users/jkoda/Research/opt/gcc/gsl
BOOST_DIR ?=
HDF5P_DIR ?=

DIR_PATH   = $(FFTW3_DIR) $(GSL_DIR) $(BOOST_DIR) $(HDF5P_DIR)

CPPFLAGS  += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS      += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

OBJS := main.o comm.o input_power.o msg.o grid.o lognormal.o
OBJS += hdf5_write.o growth.o mass_assignment.o power_spectrum.o

#
# Linking libraries
#
# LIBS += -llua -ldl
LIBS += -lboost_program_options
LIBS += -lgsl -lgslcblas -lhdf5

ifeq (,$(findstring -DDOUBLEPRECISION, $(OPT)))
  # Single precision FFTW
  FFTWSUF=f
endif
LIBS += -lfftw3$(FFTWSUF) -lfftw3$(FFTWSUF)_mpi

ifdef OPENMP
  LIBS += -lfftw3$(FFTWSUF)_omp
  #LIBS += -lfftw3$(FFTWSUF)_threads # for thread parallelization instead of omp
endif

#
# 
#
mock_lognormal_mpi: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@

# Dependences
$(OBJS): Makefile

# g++ -MM -MG *.cpp

.PHONY: clean run dependence
clean:
	rm -f $(LIB) $(OBJS)

run:
	mpirun -n 2 mock_lognormal_mpi

dependence:
	g++ -MM -MG *.cpp

