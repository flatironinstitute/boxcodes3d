for ndeg in 3 5 7; do \
   for iprec in 0 1 2 3 4; do \
     ./int2-opt-h3dtabref $iprec $ndeg
     mv fort.13 ./res/res_mac_"$ndeg"_"$iprec".txt;
   done; \
done
