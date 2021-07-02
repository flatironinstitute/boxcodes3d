EXEC = int2-voltab
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


OBJECTS = test_voltab3d.o \
    $(COM)/legetens.o \
    $(COM)/fakepolya3d.o \
    $(COM)/hkrand.o \
    $(COM)/dlaran.o \
    $(HELM)/h3danti.o \
    $(COM)/squarearbq.o \
    $(COM)/quadintrouts.o \
    $(COM)/cubeintrouts.o \
    $(COM)/voltab3d.o \
    $(COM)/qrdecomp_routs.o \
    $(COM)/loadsyms3d.o \
    $(HELM)/h3dtab.o \
    $(HELM)/h3dtab_brute.o \
    $(HELM)/helm_vol_kernels.o \

all: $(OBJECTS) 
	$(FC) $(FFLAGS)  -o $(EXEC) $(OBJECTS) $(FEND) 
	./$(EXEC) 1 3


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@

clean: 
	rm -f $(OBJECTS) $(PROJECT) fort.13
