# Makefile for fmm3dbie
# # This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create 
# the file make.inc, which overrides the defaults below (which are 
# for ubunutu linux/gcc system). 

# compiler, and linking from C, fortran
CC = gcc
CXX = g++
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy 

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp 

BOX_INSTALL_DIR=$(PREFIX)
ifeq ($(PREFIX),)
	BOX_INSTALL_DIR = ${HOME}/lib
endif

FMM_INSTALL_DIR=$(PREFIX_FMM)
ifeq ($(PREFIX_FMM),)
	FMM_INSTALL_DIR=${HOME}/lib
endif

LBLAS = -lblas -llapack

LIBS = -lm
DYLIBS = -lm
F2PYDYLIBS = -lm -lblas -llapack

LIBNAME=libboxcodes3d
DYNAMICLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a
LIMPLIB = $(DYNAMICLIB)

LFMMLINKLIB = -lfmm3d
LLINKLIB = -lboxcodes3d


# For your OS, override the above by placing make variables in make.inc
-include make.inc

# update libs and dynamic libs to include appropriate versions of
# fmm3d
#
# Note: the static library is used for DYLIBS, so that fmm3d 
# does not get bundled in with the fmm3dbie dynamic library
#
LIBS += -L$(FMM_INSTALL_DIR) $(LFMMLINKLIB) 
DYLIBS += -L$(FMM_INSTALL_DIR) $(LFMMLINKLIB)
F2PYDYLIBS += -L$(FMM_INSTALL_DIR) $(LFMMLINKLIB)

# multi-threaded libs & flags needed
ifneq ($(OMP),OFF)
  FFLAGS += $(OMPFLAGS)
  LIBS += $(OMPLIBS)
  DYLIBS += $(OMPLIBS)
  F2PYDYLIBS += $(OMPLIBS)
endif

LIBS += $(LBLAS) $(LDBLASINC)
DYLIBS += $(LBLAS) $(LDBLASINC)



# objects to compile
#
# Common objects
COM = src/Common
COMOBJS = $(COM)/cubeintrouts.o $(COM)/fakepolya3d.o \
	$(COM)/legetens.o $(COM)/loadsyms3d.o \
	$(COM)/qrdecomp_routs.o $(COM)/quadintrouts.o \
	$(COM)/rotmat_gmres.o $(COM)/squarearbq.o \
	$(COM)/tree_vol_coeffs.o $(COM)/voltab3d.o $(COM)/zerrf.o \
	$(COM)/hkrand.o $(COM)/dlaran.o


# Helmholtz wrappers
HELM = src/Helmholtz
HOBJS = $(HELM)/h3danti.o $(HELM)/h3dtab.o \
	$(HELM)/h3dtab_brute.o $(HELM)/h3dvol.o \
	$(HELM)/helm_vol_kernels.o $(HELM)/helm_volfmm3d.o \
	$(HELM)/ls_solver.o $(HELM)/helm_volfmm3d_wrap.o 

OBJS = $(COMOBJS) $(HOBJS) 

.PHONY: usage lib install test test-dyn python 

default: usage

usage:
	@echo "-------------------------------------------------------------------------"
	@echo "Makefile for fmm3dbie. Specify what to make:"
	@echo "  make install - compile and install the main library"
	@echo "  make install PREFIX=(INSTALL_DIR) - compile and install the main library at custom location given by PREFIX"
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make test - compile and run validation tests (will take around 30 secs)"
	@echo "  make test-dyn - test successful installation by validation tests linked to dynamic library (will take a couple of mins)"
	@echo "  make python - compile and test python interfaces using python"
	@echo "  make objclean - removal all object files, preserving lib & MEX"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"
	@echo ""
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=ON' for multi-threaded"
	@echo "-------------------------------------------------------------------------"



#
# implicit rules for objects (note -o ensures writes to correct dir)
#
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@



#
# build the library...
#
lib: $(STATICLIB) $(DYNAMICLIB)
ifneq ($(OMP),OFF)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif

$(STATICLIB): $(OBJS) 
	ar rcs $(STATICLIB) $(OBJS)
	mv $(STATICLIB) lib-static/

$(DYNAMICLIB): $(OBJS) 
	$(FC) -shared -fPIC $(FFLAGS) $(OBJS) -o $(DYNAMICLIB) $(DYLIBS) 
	mv $(DYNAMICLIB) lib/
	[ ! -f $(LIMPLIB) ] || mv $(LIMPLIB) lib/

install: $(STATICLIB) $(DYNAMICLIB)
	echo $(BOX_INSTALL_DIR)
	mkdir -p $(BOX_INSTALL_DIR)
	cp -f lib/$(DYNAMICLIB) $(BOX_INSTALL_DIR)/
	cp -f lib-static/$(STATICLIB) $(BOX_INSTALL_DIR)/
	[ ! -f lib/$(LIMPLIB) ] || cp lib/$(LIMPLIB) $(BOX_INSTALL_DIR)/
	@echo "Make sure to include " $(BOX_INSTALL_DIR) " in the appropriate path variable"
	@echo "    LD_LIBRARY_PATH on Linux"
	@echo "    PATH on windows"
	@echo "    DYLD_LIBRARY_PATH on Mac OSX (not needed if default installation directory is used"
	@echo " "
	@echo "In order to link against the dynamic library, use -L"$(BOX_INSTALL_DIR)  " "$(LLINKLIB) " -L"$(FMM_INSTALL_DIR)  " "$(LFMMLINKLIB)


#
# testing routines
#
test: $(STATICLIB) test/com-cube test/com-quad test/com-legediff test/com-legepq test/com-qr test/helm-legdown test/helm-legup test/helm-vol 
	cd test/Common; ./int2-cube
	cd test/Common; ./int2-quad
	cd test/Common; ./int2-legediff
	cd test/Common; ./int2-legepq
	cd test/Common; ./int2-qr
	cd test/Helmholtz; ./int2-anti-down
	cd test/Helmholtz; ./int2-anti-up
	cd test/Helmholtz; ./int2-vol
	cat print_testres.txt
	rm print_testres.txt

test-dyn: $(DYNAMICLIB) test/com-cube-dyn test/com-quad-dyn test/com-legediff-dyn test/com-legepq-dyn test/com-qr-dyn test/helm-legdown-dyn test/helm-legup-dyn test/helm-vol-dyn 
	cd test/Common; ./int2-cube
	cd test/Common; ./int2-quad
	cd test/Common; ./int2-legediff
	cd test/Common; ./int2-legepq
	cd test/Common; ./int2-qr
	cd test/Helmholtz; ./int2-anti-down
	cd test/Helmholtz; ./int2-anti-up
	cd test/Helmholtz; ./int2-vol
	cat print_testres.txt
	rm print_testres.txt

#
# Common tests static linking
#
test/com-cube: 
	$(FC) $(FFLAGS) test/Common/test_cubeintrouts.f -o test/Common/int2-cube lib-static/$(STATICLIB) $(LIBS) 

test/com-quad: 
	$(FC) $(FFLAGS) test/Common/test_quadintrouts.f -o test/Common/int2-quad lib-static/$(STATICLIB) $(LIBS) 

test/com-legediff: 
	$(FC) $(FFLAGS) test/Common/test_legediff_3d.f -o test/Common/int2-legediff lib-static/$(STATICLIB) $(LIBS) 

test/com-legepq: 
	$(FC) $(FFLAGS) test/Common/test_legetens_pqeval.f -o test/Common/int2-legepq lib-static/$(STATICLIB) $(LIBS) 

test/com-qr: 
	$(FC) $(FFLAGS) test/Common/test_dgeqp3.f90 -o test/Common/int2-qr lib-static/$(STATICLIB) $(LIBS) 

test/helm-legdown: 
	$(FC) $(FFLAGS) test/Helmholtz/test_h3danti_legedown.f -o test/Helmholtz/int2-anti-down lib-static/$(STATICLIB) $(LIBS) 

test/helm-legup: 
	$(FC) $(FFLAGS) test/Helmholtz/test_h3danti_legeup.f -o test/Helmholtz/int2-anti-up lib-static/$(STATICLIB) $(LIBS) 

test/helm-vol: 
	$(FC) $(FFLAGS) test/Helmholtz/test_helm_volfmm3d.f -o test/Helmholtz/int2-vol lib-static/$(STATICLIB) $(LIBS) 

#
# Common tests dynamic linking
#
test/com-cube-dyn: 
	$(FC) $(FFLAGS) test/Common/test_cubeintrouts.f -o test/Common/int2-cube -L$(FMM_INSTALL_DIR) -L$(BOX_INSTALL_DIR) $(LFMMLINKLIB) $(LLINKLIB) $(LIBS)

test/com-quad-dyn: 
	$(FC) $(FFLAGS) test/Common/test_quadintrouts.f -o test/Common/int2-quad -L$(FMM_INSTALL_DIR) -L$(BOX_INSTALL_DIR) $(LFMMLINKLIB) $(LLINKLIB) $(LIBS)

test/com-legediff-dyn: 
	$(FC) $(FFLAGS) test/Common/test_legediff_3d.f -o test/Common/int2-legediff -L$(FMM_INSTALL_DIR) -L$(BOX_INSTALL_DIR) $(LFMMLINKLIB) $(LLINKLIB) $(LIBS)

test/com-legepq-dyn: 
	$(FC) $(FFLAGS) test/Common/test_legetens_pqeval.f -o test/Common/int2-legepq -L$(FMM_INSTALL_DIR) -L$(BOX_INSTALL_DIR) $(LFMMLINKLIB) $(LLINKLIB) $(LIBS)

test/com-qr-dyn: 
	$(FC) $(FFLAGS) test/Common/test_dgeqp3.f90 -o test/Common/int2-qr -L$(FMM_INSTALL_DIR) -L$(BOX_INSTALL_DIR) $(LFMMLINKLIB) $(LLINKLIB) $(LIBS)

test/helm-legdown-dyn: 
	$(FC) $(FFLAGS) test/Helmholtz/test_h3danti_legedown.f -o test/Helmholtz/int2-anti-down -L$(FMM_INSTALL_DIR) -L$(BOX_INSTALL_DIR) $(LFMMLINKLIB) $(LLINKLIB) $(LIBS)

test/helm-legup-dyn: 
	$(FC) $(FFLAGS) test/Helmholtz/test_h3danti_legeup.f -o test/Helmholtz/int2-anti-up -L$(FMM_INSTALL_DIR) -L$(BOX_INSTALL_DIR) $(LFMMLINKLIB) $(LLINKLIB) $(LIBS)

test/helm-vol-dyn: 
	$(FC) $(FFLAGS) test/Helmholtz/test_helm_volfmm3d.f -o test/Helmholtz/int2-vol -L$(FMM_INSTALL_DIR) -L$(BOX_INSTALL_DIR) $(LFMMLINKLIB) $(LLINKLIB) $(LIBS)


#
# build the python bindings/interface
#
python: $(STATICLIB)
	cd python && export FMMBIE_LIBS='$(LIBS)' && pip install -e . 


#
# build python gmsh to go3
#
python-gmsh: $(DYNAMICLIB)
	f2py $(F2PYDYLIBS) -lfmm3dbie -c src/tria_routs/koornexps.f90 -m kexp
	f2py $(F2PYDYLIBS) -lfmm3dbie -c src/surface_routs/surf_routs.f90 -m srout
	mv kexp*.so python/
	mv srout*.so python/
	rm -rf *.dSYM

#
# housekeeping routines
#
clean: objclean
	rm -f lib-static/*.a lib/*.so
	rm -f test/Common/int2*
	rm -f test/Helmholtz/int2*
	rm -f python/*.so
	rm -rf python/build
	rm -rf python/fmm3dpy.egg-info

objclean: 
	rm -f $(OBJS) $(TOBJS)
	rm -f test/Helmholtz/*.o test/Common/*.o 
