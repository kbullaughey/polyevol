#!/bin/bash

#$ -cwd
#$ -l h_vmem=1g
#$ -o grid/equilibrium.$JOB_ID.out
#$ -e grid/equilibrium.$JOB_ID.err

if [ "x" == "x$1" ] ; then
  echo "Must specify effects" >&2
  exit 1
fi
effects=$*

seed=`$HOME/bin/rand.pl`
file_suffix=`echo "$effects" | sed 's/--/-/g' | tr '= ' '_-' | sed 's/^-//'`
outfile="out/$file_suffix/equilibrium-$file_suffix.$JOB_ID.out"

if [ ! -d out/$file_suffix ] ; then
  mkdir out/$file_suffix
fi

$HOME/bin/quant -m infinite --freqs=even --seed=$seed -s 0.1 -N 1000 -u 1e-03 \
    $effects --opts=0 --loci=20 --burnin=5000 --times=500000 \
 | awk '
    BEGIN {FS="[: ]";} 
    /^site: id:/ {effects[$5] = $0; next; } 
    $4 == "freqs" && $3 % 100 == 0 { 
      for (i=6; i < NF; i+=2) { 
        if ($i in effects) { 
          print effects[$i]; 
          delete effects[$i]; 
        } 
      } 
    } 
    $3 % 100 == 0 { print $0; }' \
 > $outfile

R --vanilla --args --file=$outfile < save_frequency_distribution.r

# END
