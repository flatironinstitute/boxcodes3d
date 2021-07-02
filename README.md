# Boxcodes

## Fortran routines for evaluating volume potentials for Laplace and Helmholtz equation

This repository contains fmm acclerated codes for evaluating volume potentials for Helmholtz equations
and iterative solvers for the solution of volumetric scattering problems like the Lippmann Schwinger
equation.

The repository has an external dependency - [FMM3D](https://fmm3d.readthedocs.io/en/latest)

For easy installation, choose an appropriate make.inc file based on your operating
system and compiler,  and run ``make install'' in the FMM3D directory followed by ``make install''
in this repository.

To verify successful installation run ``make test" or ``make test-dyn". 
Note that linking to the shared object in ``make test-dyn" will require
the LD_LIBRARY_PATH in linux to be appropriately set.

To see examples of using the volume potential see the examples folder.
The example folder demos both the evaluation of the volume potential
and a demo for the lippman schwinger solver.

The near quadrature for the volume potential is computed using an
appropriate inverse PDE solver coupled with adaptive integration
on the surface of the unit cube. Since the kernels are translationally
invariant, this computation needs to be done only once per level.
The adaptive integration routines have been optimized for improved
performance for this task. The far field is accelerated using fast
multipole methods.


##Upcoming

- Support for Maxwell and Stokes volume potentials
- Python interfaces
- Julia interfaces
- Matlab interfaces  
