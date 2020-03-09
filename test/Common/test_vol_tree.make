
EXEC = vol_tree

#HOST = osx
HOST=linux-gfortran
#HOST=linux-ifort
#HOST=linux-gfortran-prof
HOST=linux-gfortran-openmp

ifeq ($(HOST),osx)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -c -w
FLINK = gfortran -w -o $(EXEC)
FEND = -framework accelerate
endif

ifeq ($(HOST),linux-gfortran)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math -c -w  
FLINK = gfortran -w -o $(EXEC) 
FEND = -lopenblas 
endif

ifeq ($(HOST),linux-gfortran-prof)
FC = gfortran
FFLAGS = -O3 -march=native -pg -g -funroll-loops -ftree-vectorize -ffast-math -c -w  
FLINK = gfortran -w -o $(EXEC) -pg
FEND = -lblas -llapack
endif

ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran
FFLAGS = -O3 -march=native -funroll-loops -ftree-vectorize -ffast-math -c -w  --openmp
FLINK = gfortran -w --openmp -o $(EXEC) 
FEND = -lopenblas -llapack
endif

ifeq ($(HOST),linux-ifort)
FC = ifort
FFLAGS = -O1 -g -c -w -xW -qopenmp 
FLINK = ifort -w -qopenmp -o $(EXEC)
WITH_SECOND = 1
endif


SRC = ../../src
UTILS_DIR = ../../../utils


.PHONY: all clean list

SOURCES =  test_vol_tree.f \
  $(SRC)/Common/prini_new.f \
  $(UTILS_DIR)/legeexps.f \
  $(SRC)/Common/tree_vol.f \
  $(SRC)/Common/legetens.f \
  $(SRC)/Common/voltab3d.f \
  $(UTILS_DIR)/hkrand.f \
  $(UTILS_DIR)/dlaran.f \

ifeq ($(WITH_SECOND),1)
SOURCES += $(SRC)/second-r8.f
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
	./$(EXEC) 2 

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



