#!/bin/bash

#$ -cwd
#$ -l h_vmem=2g
#$ -o grid/equilibrium.$JOB_ID.out
#$ -e grid/equilibrium.$JOB_ID.err

if [ "x" == "x$1" ] ; then
  echo "Must specify params" >&2
  exit 1
fi
params=$*
analysis="phenotype-N_100-s_0.005"

seed=`$HOME/bin/rand.pl`
if [ -z $JOB_ID ] ; then
  JOB_ID=$seed
fi

file_suffix=`echo "$params" | tr '= ' '_-' | sed 's/--/-/g' | sed 's/^-//'`
basedir="out/$analysis"
outfile="$basedir/$analysis-$file_suffix.$JOB_ID.out"

if [ ! -d $basedir ] ; then
  mkdir -p $basedir
fi

$HOME/bin/quant -m infinite --freqs=even -s 0.005 -N 100 --opts=0 --loci=0 --burnin=2000 \
    --seed=$seed --times=1000 --disable-all-stats --enable-stat=phenotype $params \
  > $outfile

# END
