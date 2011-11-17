#!/bin/bash

seed=1
prog="../quant"
cmd="cat qin20 | $prog -m infinite -N 1000 -u 0.001 -s 0.12 --seed=$seed --burnin=5000 --opt=0,2.0 --effects=1 --eprob=1 --times=1000,2000 --loci=20 > run.out"

echo "beginning run"
echo $cmd
eval $cmd
echo "plotting"
if R --vanilla < plot.r >/dev/null; then 
  echo "done"
else 
  echo "failed"
fi
