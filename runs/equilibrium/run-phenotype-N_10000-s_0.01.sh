#!/bin/bash

#$ -cwd
#$ -l h_vmem=4g
#$ -o grid/equilibrium.$JOB_ID.out
#$ -e grid/equilibrium.$JOB_ID.err

if [ "x" == "x$1" ] ; then
  echo "Must specify params" >&2
  exit 1
fi
params=$*
analysis="phenotype-N_10000-s_0.01"

seed=`$HOME/bin/rand.pl`
if [ -z $JOB_ID ] ; then
  JOB_ID=$seed
fi

file_suffix=`echo "$params" | tr '= ' '_-' | sed 's/--/-/g' | sed 's/^-//'`
basedir="out/$analysis"
outfile="$basedir/$analysis-$file_suffix.$JOB_ID.out.gz"

if [ ! -d $basedir ] ; then
  mkdir -p $basedir
fi

$HOME/bin/quant -m infinite --freqs=even -s 0.01 -N 10000 --opts=0 --loci=0 --burnin=4000 \
    --times=10000 --disable-all-stats --enable-stat=phenotype $params \
  | gzip > $outfile

# END
