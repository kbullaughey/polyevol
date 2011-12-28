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

seed=`$HOME/bin/rand.pl`
file_suffix=`echo "$params" | tr '= ' '_-' | sed 's/--/-/g' | sed 's/^-//'`
outfile="out/sojourn/$file_suffix/equilibrium-$file_suffix.$JOB_ID.out.gz"

if [ ! -d out/sojourn/$file_suffix ] ; then
  mkdir -p out/sojourn/$file_suffix
fi

$HOME/bin/quant -m infinite --freqs=even --seed=$seed -s 0.1 -N 1000 \
    $params --opts=0 --loci=20 --burnin=5000 --times=100000 \
 | gzip > $outfile

zcat $outfile \
 | awk '/end burnin/{on=1} { if(on) { print $0; } }' \
 | grep 'absorption' \
 | sed 's/^.*sojourn: \([0-9]*\) effect: \(.*\)$/\1 \2/' \
 | gzip \
 > `echo $outfile | sed 's/.out.gz/.sojourns.gz/'`

#R --vanilla --args --file=$outfile < save_frequency_distribution.r

# END
