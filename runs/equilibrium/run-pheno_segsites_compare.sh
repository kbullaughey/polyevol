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
analysis="pheno_segsites_compare"

seed=`$HOME/bin/rand.pl`
if [ -z $JOB_ID ] ; then
  JOB_ID=$seed
fi

file_suffix=`echo "$params" | tr '= ' '_-' | sed 's/--/-/g' | sed 's/^-//'`
basedir="out/$analysis/$file_suffix"
outfile="$basedir/$analysis-$file_suffix.$JOB_ID.out"

if [ ! -d $basedir ] ; then
  mkdir -p $basedir
fi

$HOME/bin/quant -m infinite --freqs=even -s 0.1 -N 1000 --opts=0 --loci=0 --burnin=2000 \
    --effects=1,0 --eprobs=0.5,0.5 --times=10000 --disable-all-stats \
    --enable-stat=phenotype --enable-stat=segsites $params \
  | awk '
    NR <= 2 {
      print $0;
    } NR % 100 == 0 {
      print NR > "/dev/stderr"
    } /segsites:/ {
      for (i=4; i<NF; i++) { 
        split($i, pair, ","); 
        if (pair[1] ~ /^-*0$/) { neutral += pair[2]; } 
      } 
      neutral_count++; 
    } /pheno:/ {
      phenovarsum += $5;
      phenovarcount ++;
    } END {
      print "neutral:", neutral, neutral_count;
      print "phenovar:", phenovarsum, phenovarcount;
    }
  '> $outfile
# END
