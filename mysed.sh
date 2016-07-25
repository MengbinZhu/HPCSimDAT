#!/bin/bash
 var=qgcont.f90
  echo "Now the f90 is:      ${var}"
  for ((xx=1;xx<=50;xx=xx+10)) ; do
    cp ${var} a${xx}_${var}
    echo "Now the j is: a${xx}_${var}"
    sed -i "s#isample=1,nsample#isample=${xx},${xx}+9#g" a${xx}_$var
    #ifort -o $xx.out -I /disk1/soft/netcdf-3.6.0-x64/include -L /disk1/soft/netcdf-3.6.0-x64/lib -lnetcdf $var
    ifort -o $xx.out a${xx}_$var
#cp sub.sh sub$xx.sh
#sed -i "s#a.out#${xx}.out#g" sub$xx.sh
    nohup ./$xx.out >fj_data.$xx &
    echo "qsub ok"
  done
