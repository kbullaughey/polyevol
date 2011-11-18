#!/bin/bash

if [ ! -d out ] ; then
  mkdir out
fi

s=0.1

for i in `seq 1 4`; do
  echo "iteration $i"
  quant -m infinite --freqs=even --seed=`rand.pl` -s $s -N 1000 -u 1e-03 \
    --effects=1 --opts=0,1 --loci=20 --burnin=2000 --times=100,10900 \
    > out/shift_opt_0_to_1-times_100_10900-s_$s.$i.out
done

# END
