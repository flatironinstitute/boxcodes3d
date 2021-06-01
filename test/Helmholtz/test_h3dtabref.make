
EXEC = int2-opt-h3dtabref 

HOST = osx
#HOST=linux-gfortran
#HOST=linux-gfortran-openmp
#HOST=linux-gfortran-debug

ifeq ($(HOST),osx)
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -funroll-loops -c -w -fopenmp -std=legacy
FLINK = gfortran -w -fopenmp -std=legacy -o $(EXEC)
FEND = -lopenblas ${LDFLAGS}
endif

ifeq ($(HOST),linux-gfortran)
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math -c -w  
FLINK = gfortran -w -o $(EXEC) 
FEND = -lopenblas 
endif

ifeq ($(HOST),linux-gfortran-debug)
FC = gfortran
FFLAGS = -g -c -w  
FLINK = gfortran -w -g -o $(EXEC) 
FEND = -lblas -llapack
endif

ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -fopenmp -funroll-loops -c -w  
FLINK = gfortran -w -fopenmp -o $(EXEC) 
FEND = -lopenblas
endif

ifeq ($(HOST),linux-ifort)
FC = ifort
FFLAGS = -O3 -c -w -xW 
FLINK = ifort -w -mkl -o $(EXEC)
WITH_SECOND = 1
endif


SRC = ../../src


.PHONY: all clean list

SOURCES =  test_h3dtabref.f \
  $(SRC)/Common/prini_new.f \
  $(SRC)/Common/legeexps.f \
  $(SRC)/Common/legetens.f \
  $(SRC)/Common/squarearbq.f \
  $(SRC)/Common/hkrand.f \
  $(SRC)/Common/dlaran.f \
  $(SRC)/Common/sphere_pol_routs.f \
  $(SRC)/Common/ncleastsq.f \
  $(SRC)/Common/svdpivot.f \
  $(SRC)/Common/csvdpiv.f \
  $(SRC)/Common/qleigen_trid.f \
  $(SRC)/Common/yrecursion.f \
  $(SRC)/Common/quadintrouts2.f \
  $(SRC)/Common/quadintrouts.f \
  $(SRC)/Common/voltab3d.f \
  $(SRC)/Helmholtz/lommel.f \
  $(SRC)/Helmholtz/h3dtab.f \
  $(SRC)/Helmholtz/h3dtab_brute.f \
  $(SRC)/Common/cubeintrouts2.f


ifeq ($(WITH_SECOND),1)
SOURCES += $(HELLSKITCHEN)/Common/second-r8.f
endif

OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))

#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f
	$(FC) $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

all: $(OBJECTS)
	rm -f $(EXEC)
	$(FLINK) $(OBJECTS) $(FEND)
	./$(EXEC) 1 3 

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



