#!/bin/bash

#$ -cwd
#$ -l h_vmem=1g
#$ -o grid/equilibrium.$JOB_ID.out
#$ -e grid/equilibrium.$JOB_ID.err

if [ "x" == "x$1" ] ; then
  echo "Must specify params" >&2
  exit 1
fi
params=$*
analysis="visits-2012_01_31"

seed=`$HOME/bin/rand.pl`
file_suffix=`echo "$params" | tr '= ' '_-' | sed 's/--/-/g' | sed 's/^-//'`
basedir="out/$analysis/$file_suffix"
outfile="$basedir/$analysis-$file_suffix.$JOB_ID.out.gz"

if [ ! -d $basedir ] ; then
  mkdir -p $basedir
fi

$HOME/bin/quant -m infinite --freqs=even -s 0.1 -N 1000 --opts=0 --loci=0 --burnin=10000 --seed=$seed \
    --times=1000000 --disable-all-stats --enable-stat=visits $params \
  | gzip > $outfile

# END
