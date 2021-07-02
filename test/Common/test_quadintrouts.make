EXEC = int2-quad
#HOST = gcc
HOST = gcc-openmp
#HOST = intel
#HOST = intel-ompenmp

#
# For linux systems, it is assumed that the environment
# variable LD_LIBRARY_PATH contains the locations to libfmm3d.so
# for Macosx, the .so file also need to be copied over /usr/local/lib


FMM_INSTALL_DIR = $(PREFIX_FMM)
LBLAS = $(PREFIX_BLAS)
LBLASINC = $(PREFIX_BLAS_INC)
LFMMLINKLIB = -lfmm3d

ifneq ($(OS),Windows_NT) 
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        ifeq ($(PREFIX_FMM),)
            FMM_INSTALL_DIR=/usr/local/lib
        endif

        ifeq ($(PREFIX_BLAS),)
            LBLAS = -framework accelerate
        endif
    endif
    ifeq ($(UNAME_S),Linux)
        ifeq ($(PREFIX_FMM),)
            FMM_INSTALL_DIR=${HOME}/lib
        endif
        ifeq ($(PREFIX_BLAS),)
            LBLAS = -lblas -llapack
        endif
    endif
endif


ifeq ($(HOST),gcc)
    FC=gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -std=legacy 
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp -std=legacy 
endif

ifeq ($(HOST),intel)
    FC=ifort 
    FFLAGS= -O3 -fPIC -march=native
endif

ifeq ($(HOST),intel-openmp)
    FC = ifort 
    FFLAGS= -O3 -fPIC -march=native -qopenmp
endif

FEND = -L${FMM_INSTALL_DIR} $(LFMMLINKLIB) $(LBLAS) $(LDBLASINC)



# Test objects
#
COM = ../../src/Common
HELM = ../../src/Helmholtz

.PHONY: all clean

default: all


OBJECTS = test_quadintrouts.o \
    $(COM)/legetens.o \
    $(COM)/squarearbq.o \
    $(COM)/hkrand.o \
    $(COM)/dlaran.o \
    $(COM)/quadintrouts.o \

all: $(OBJECTS) 
	$(FC) $(FFLAGS)  -o $(EXEC) $(OBJECTS) $(FEND) 
	./$(EXEC)


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

clean: 
	rm -f $(OBJECTS) $(PROJECT) fort.13
