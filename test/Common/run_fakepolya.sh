#!/bin/bash
#
# test runner

mkdir -p output
rm -f int2-fakepolya
make -f test_fakepolya.make
rm -f output/fakepolya*.txt

for i in {4..16..2}
do
    for j in {1,1.25,1.5,1.75,2}
    do
	echo "$i" > tmp.txt
	
	echo "$j" >> tmp.txt
	./int2-fakepolya < tmp.txt
    done
done	      

rm -f tmp.txt
