#!/bin/bash

# Here we consider the number of segregating sites for various mutation rates 
#   at equilibrium with no environmental change

for i in `seq 1 64`; do 
#  qsub run-segregating_sites.sh --effects=0 -u ${i}e-03
  qsub run-segregating_sites.sh --effects=0.5 -u ${i}e-03
  qsub run-segregating_sites.sh --effects=1 -u ${i}e-03
done

# END
