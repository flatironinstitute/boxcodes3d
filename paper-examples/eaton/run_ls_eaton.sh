#! /bin/bash

make -f test_ls_eaton.make

echo "./" > infile
echo "$(date --utc +%Y%m%d_%H%M%S)" >> infile
echo "1" >> infile
echo "2" >> infile
echo "1e-4" >> infile
echo "1" >> infile

./int2-fmm < infile
