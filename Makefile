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

DIR_PATH   = $(FFTW3_DIR) $(GSL_DIR)

CPPFLAGS  += $(foreach dir, $(DIR_PATH), -I$(dir)/include)
LIBS      += $(foreach dir, $(DIR_PATH), -L$(dir)/lib)

OBJS := main.o #comm.o msg.o config.o fft.o mem.o particle.o util.o
OBJS += #power.o cosmology.o lpt.o

#
# Linking libraries
#
# LIBS += -llua -ldl 
LIBS += -lgsl -lgslcblas

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
	$(CXX) $(OBJS) -o $@

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

