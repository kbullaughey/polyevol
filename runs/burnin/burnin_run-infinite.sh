#!/bin/bash

#$ -cwd
#$ -l h_vmem=800m
#$ -o grid/burnin.$JOB_ID.out
#$ -e grid/burnin.$JOB_ID.err

quant="$HOME/bin/quant"
seed=`$HOME/bin/rand.pl`

if [ "x" == "x$1" ] ; then
  echo "must specify s" >&2
  exit 1
fi
s=$1

if [ "x" == "x$2" ] ; then
  echo "must specify reps" >&2
  exit 1
fi
reps=$2

for i in `seq 1 $reps`; do 
  if [ ! -d out-infinite/$s/$i ] ; then
    mkdir -p out-infinite/$s/$i
  fi
  $quant -m infinite -s $s --freqs=even -N 1000 --times=1000 --seed=$seed --opt=1 --effects=1 \
    --loci=25 --burnin=0 > out-infinite/$s/$i/burnin.$JOB_ID.out
done

# END
