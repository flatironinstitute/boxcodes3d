for icase in 3; do \
   for norder in 4 6 8 12; do \
      for iprec in 0 1 2 3 4; do \
         ./int2-volvspt $iprec $norder $icase
         mv fort.13 ./timing-res/res_linux_"$norder"_"$iprec"_"$icase".txt;
      done; \
   done; \
done
