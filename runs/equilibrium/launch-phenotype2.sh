#!/bin/bash

for i in `cat mu_0.001-0.15_by_0.0051`; do 
  qsub run-phenotype.sh -u $i --effects=1
done
