#!/bin/bash

for i in `cat 20_mutation_rates`; do 
  qsub run-phenotype.sh -u $i --effects=1
done
