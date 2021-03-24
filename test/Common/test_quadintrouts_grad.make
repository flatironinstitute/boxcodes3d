
EXEC = triaintrouts

#HOST = osx
HOST=linux-gfortran

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

SOURCES =  test_quadintrouts_grad.f \
  $(SRC)/Common/prini_new.f \
  $(SRC)/Common/legeexps.f \
  $(SRC)/Common/legetens.f \
  $(SRC)/Common/squarearbq.f \
  $(SRC)/Common/hkrand.f \
  $(SRC)/Common/dlaran.f \
  $(SRC)/Common/quadintrouts.f \

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



