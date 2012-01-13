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
analysis="segregating_sites"

seed=`$HOME/bin/rand.pl`
file_suffix=`echo "$params" | tr '= ' '_-' | sed 's/--/-/g' | sed 's/^-//'`
outfile="out/$analysis/$file_suffix/equilibrium-$file_suffix.$JOB_ID.out.gz"
summaryfile="out/$analysis/$file_suffix/segsite_counts-$file_suffix.$JOB_ID.out.gz"

if [ ! -d out/$analysis/$file_suffix ] ; then
  mkdir -p out/$analysis/$file_suffix
fi

$HOME/bin/quant -m infinite --freqs=even --seed=$seed -s 0.1 -N 1000 \
    $params --opts=0 --loci=20 --burnin=10000 --times=100000 --disable-all-stats --enable-stat=frequencies \
 | gzip > $outfile

zcat $outfile | awk 'NR <= 2 { print $0; } /freqs:/ {print NF-3;}' | gzip > $summaryfile
rm $outfile

# END
