#!/bin/bash

#$ -cwd
#$ -l h_vmem=1g
#$ -o grid/equilibrium.$JOB_ID.out
#$ -e grid/equilibrium.$JOB_ID.err

if [ "x" == "x$1" ] ; then
  echo "Must specify N" >&2
  exit 1
fi
N=$1
if [ "x" == "x$2" ] ; then
  echo "Must specify Nu" >&2
  exit 1
fi
Nu=$2
if [ "x" == "x$3" ] ; then
  echo "Must specify Nsa^2" >&2
  exit 1
fi
Nsa2=$3

# compute the individual parameters from the composites, I assume an effect size 1
u=`python -c "print (0.0+$Nu)/$N"`
s=`python -c "print (0.0+$Nsa2)/$N"`

analysis='fix_composite'
seed=`$HOME/bin/rand.pl`
file_suffix="N_$N-Nu_$Nu-Nsa2_$Nsa2"
outfile="out/$analysis/$file_suffix/equilibrium-$file_suffix.$JOB_ID.out.gz"

if [ ! -d out/$analysis/$file_suffix ] ; then
  mkdir -p out/$analysis/$file_suffix
fi

$HOME/bin/quant -m infinite --freqs=even --seed=$seed --effects=1 \
    -N $N -u $u -s $s --opts=0 --loci=20 --burnin=5000 --times=100000 \
 | gzip > $outfile

zcat $outfile \
 | awk '/end burnin/{on=1} { if(on) { print $0; } }' \
 | grep 'absorption' \
 | sed 's/^.*sojourn: \([0-9]*\) effect: \(.*\)$/\1 \2/' \
 | gzip \
 > `echo $outfile | sed 's/.out.gz/.sojourns.gz/'`

ssitefile=`echo $outfile | sed 's/out.gz$/segsites/'`
zcat $outfile | grep 'freqs:' | awk '{print NF-3}' | gzip > $ssitefile.gz


# END
