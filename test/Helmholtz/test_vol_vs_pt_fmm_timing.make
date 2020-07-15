
EXEC = int2-fmm

HOST = osx
HOST=linux-gfortran

FC = gfortran
FFLAGS = -fPIC -O3 -lstdc++ -march=native -fopenmp -funroll-loops
CXXFLAGS = -std=c++11 -lstdc++ -DSCTL_PROFILE=-1 -fPIC -O3 -march=native -fopenmp -funroll-loops -I../../../FMM3D/vec-kernels/include 
FLINK = gfortran -w -fopenmp -o $(EXEC)
FEND = -lopenblas ${LDFLAGS} -lstdc++


SRC = ../../src
FMM3D = ../../../FMM3D/src
UTILS_DIR = ../../../utils


.PHONY: all clean list

OBJECTS =  test_vol_vs_pt_fmm_timing.o \
  $(SRC)/Common/prini_new.o \
  $(UTILS_DIR)/legeexps.o \
  $(SRC)/Common/tree_vol.o \
  $(SRC)/Common/legetens.o \
  $(SRC)/Common/voltab3d.o \
  $(SRC)/Helmholtz/h3dvol.o \
  $(SRC)/Helmholtz/h3dtab.o \
  $(SRC)/Helmholtz/lommel.o \
  $(SRC)/Helmholtz/helm_volfmm3d.o \
  $(SRC)/Common/sphere_pol_routs.o \
  $(SRC)/Common/ncleastsq.o \
  $(SRC)/Common/svdpivot.o \
  $(SRC)/Common/csvdpiv.o \
  $(SRC)/Common/qleigen_trid.o \
  $(SRC)/Common/yrecursion.o \
  $(SRC)/Common/quadintrouts.o \
  $(SRC)/Common/quadintrouts2.o \
  $(SRC)/Common/loadsyms3d.o \
  $(SRC)/Common/squarearbq.o \
  $(SRC)/Common/zerrf.o \
  $(FMM3D)/Helmholtz/hfmm3d.o \
  $(FMM3D)/Helmholtz/hfmm3dwrap.o \
  $(FMM3D)/Helmholtz/hpwrouts.o \
  $(FMM3D)/Helmholtz/h3dtrans.o \
  $(FMM3D)/Helmholtz/h3dterms.o \
  $(FMM3D)/Helmholtz/helmrouts3d.o \
  $(FMM3D)/Helmholtz/projections.o \
  $(FMM3D)/Common/rotviarecur.o \
  $(FMM3D)/Common/rotproj.o \
  $(FMM3D)/Common/besseljs3d.o \
  $(FMM3D)/Common/tree_lr_3d.o \
  $(FMM3D)/Common/fmmcommon.o \
  $(FMM3D)/Common/rotgen.o \
  $(FMM3D)/Common/dfft.o \
  $(FMM3D)/Helmholtz/h3dcommon.o \
  $(FMM3D)/Helmholtz/hwts3e.o \
  $(FMM3D)/Helmholtz/hnumphys.o \
  $(FMM3D)/Helmholtz/hnumfour.o \
  $(UTILS_DIR)/hkrand.o \
  $(UTILS_DIR)/dlaran.o \
  $(SRC)/Common/aquad.o \
  $(SRC)/Common/cerf.o \


ifeq ($(WITH_SECOND),1)
OBJECTS += $(SRC)/second-r8.o
endif

ifeq ($(FAST_KER), ON)
  OBJECTS += $(FMM3D)/Helmholtz/helmkernels_fast.o
  OBJECTS += $(FMM3D)/Helmholtz/hndiv_fast.o
  OBJECTS += $(FMM3D)/../vec-kernels/src/libkernels.o
endif

ifneq ($(FAST_KER), ON)
  OBJECTS += $(FMM3D)/Helmholtz/helmkernels.o
  OBJECTS += $(FMM3D)/Helmholtz/hndiv.o
endif

#
# use only the file part of the filename, then manually specify
# the build location
#

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@
%.o : %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@
%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

all: $(OBJECTS)
	rm -f $(EXEC)
	$(FLINK) $(OBJECTS) $(FEND)
	./$(EXEC) 1 4 3 

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



