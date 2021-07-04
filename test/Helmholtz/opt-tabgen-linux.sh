for ndeg in 3 5 7 11; do \
   for iprec in 0 1 2 3 4; do \
     ./int2-tabref $iprec $ndeg
     mv fort.13 ./res/res_linux_"$ndeg"_"$iprec".txt;
   done; \
done
