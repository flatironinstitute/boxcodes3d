#!/bin/sh -v

gfortran -c -g ../prini_new.f ../voltab3d.f
gfortran -o int2 -g genloadsyms3d.f prini_new.o voltab3d.o

./int2
