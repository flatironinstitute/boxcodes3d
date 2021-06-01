
EXEC = int2h3danti_form

#HOST = osx
HOST=linux-gfortran
HOST=linux-gfortran-openmp
#HOST=linux-gfortran-debug

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
FEND = -lblas -llapack
endif

ifeq ($(HOST),linux-gfortran-debug)
FC = gfortran
FFLAGS = -g -c -w  
FLINK = gfortran -w -g -o $(EXEC) 
FEND = -lblas -llapack
endif

ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran
FFLAGS = -O3 -march=native --openmp -funroll-loops -ftree-vectorize -ffast-math -c -w  
FLINK = gfortran -w --openmp -o $(EXEC) 
FEND = -lblas -llapack
endif

ifeq ($(HOST),linux-ifort)
FC = ifort
FFLAGS = -O3 -c -w -xW 
FLINK = ifort -w -mkl -o $(EXEC)
WITH_SECOND = 1
endif


SRC = ../../src
UTILS_DIR = ../../utils
XTRI_DIR = ../../xtri/src
YTRI_DIR = ../../ytri/src


.PHONY: all clean list

SOURCES =  test_h3danti_form.f \
  $(SRC)/Common/prini_new.f \
  $(SRC)/Common/legeexps.f \
  $(SRC)/Common/legetens.f \
  $(SRC)/Helmholtz/h3danti.f 


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
	./$(EXEC) 2 

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)


