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
analysis="phenotype_effect_grid"

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

$HOME/bin/quant -m infinite --freqs=even -s 0.1 -N 1000 --opts=0 --loci=0 --burnin=10000 \
    --times=40000 --disable-all-stats --enable-stat=phenotype $params \
  | awk '
    /^cmd|params/ {
      print $0; 
    } NR % 100 == 0 {
      print phenovarcount > "/dev/stderr"
    } /pheno:/ {
      phenovarsum += $5;
      phenovarcount++;
    } END {
      print "phenovarsum:", phenovarsum, "count:", phenovarcount;
    }' \
  > $outfile

# END
