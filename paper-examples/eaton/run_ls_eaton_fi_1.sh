#! /bin/bash

echo "/mnt/ceph/users/taskham/" > infile
echo "$(date --utc +%Y%m%d_%H%M%S)" >> infile
echo "1" >> infile
echo "3" >> infile
echo "1e-5" >> infile
echo "1" >> infile

./int2-fmm < infile
