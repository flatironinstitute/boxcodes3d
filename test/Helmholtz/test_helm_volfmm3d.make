
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
FFLAGS = -fPIC -O3 -march=native -funroll-loops -c -w  
FLINK = gfortran -w -o $(EXEC) 
FEND = -lopenblas -L/usr/local/opt/openblas/lib 
endif

ifeq ($(HOST),linux-gfortran-prof)
FC = gfortran
FFLAGS = -O3 -march=native -pg -g -funroll-loops -ftree-vectorize -ffast-math -c -w  
FLINK = gfortran -w -o $(EXEC) -pg
FEND = -lblas -llapack
endif

ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -c --openmp
FLINK = gfortran -w --openmp -o $(EXEC) 
FEND = -lopenblas -L/usr/local/opt/openblas/lib 
endif

ifeq ($(HOST),linux-ifort)
FC = ifort
FFLAGS = -O1 -g -c -w -xW -qopenmp 
FLINK = ifort -w -qopenmp -o $(EXEC)
WITH_SECOND = 1
endif


SRC = ../../src
FMM3D = ../../../FMM3D/src
UTILS_DIR = ../../../utils


.PHONY: all clean list

SOURCES =  test_helm_volfmm3d.f \
  $(SRC)/Common/prini_new.f \
  $(UTILS_DIR)/legeexps.f \
  $(SRC)/Common/tree_vol.f \
  $(SRC)/Common/legetens.f \
  $(SRC)/Common/voltab3d.f \
  $(SRC)/Helmholtz/h3dvol.f \
  $(SRC)/Helmholtz/h3dtab.f \
  $(SRC)/Helmholtz/lommel.f \
  $(SRC)/Helmholtz/helm_volfmm3d.f \
  $(SRC)/Common/sphere_pol_routs.f \
  $(SRC)/Common/ncleastsq.f \
  $(SRC)/Common/svdpivot.f \
  $(SRC)/Common/csvdpiv.f \
  $(SRC)/Common/qleigen_trid.f \
  $(SRC)/Common/yrecursion.f \
  $(SRC)/Common/quadintrouts.f \
  $(SRC)/Common/loadsyms3d.f \
  $(SRC)/Common/squarearbq.f \
  $(FMM3D)/Helmholtz/hpwrouts.f \
  $(FMM3D)/Helmholtz/h3dtrans.f \
  $(FMM3D)/Helmholtz/h3dterms.f \
  $(FMM3D)/Helmholtz/helmrouts3d.f \
  $(FMM3D)/Helmholtz/projections.f \
  $(FMM3D)/Common/rotviarecur.f \
  $(FMM3D)/Common/rotproj.f \
  $(FMM3D)/Common/besseljs3d.f \
  $(FMM3D)/Common/fmmcommon.f \
  $(FMM3D)/Common/rotgen.f \
  $(FMM3D)/Common/dfft.f \
  $(FMM3D)/Helmholtz/h3dcommon.f \
  $(FMM3D)/Helmholtz/hwts3e.f \
  $(FMM3D)/Helmholtz/hnumphys.f \
  $(FMM3D)/Helmholtz/hnumfour.f \
  $(UTILS_DIR)/hkrand.f \
  $(UTILS_DIR)/dlaran.f \
  $(SRC)/Common/aquad.f \
  $(SRC)/Common/cerf.f90 \

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



