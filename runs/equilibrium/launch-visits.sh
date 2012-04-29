#!/bin/bash

for i in `seq 1 10`; do qsub run-visits.sh -u 1e-03 --effects=1; done
for i in `seq 1 10`; do qsub run-visits.sh -u 4e-03 --effects=1; done
for i in `seq 1 10`; do qsub run-visits.sh -u 1.6e-02 --effects=1; done
for i in `seq 1 10`; do qsub run-visits.sh -u 6.4e-02 --effects=1; done
for i in `seq 1 10`; do qsub run-visits.sh -u 3.2e-02 --effects=1; done
for i in `seq 1 10`; do qsub run-visits.sh -u 8e-03 --effects=1; done
for i in `seq 1 10`; do qsub run-visits.sh -u 1.28e-01 --effects=1; done
for i in `seq 1 10`; do qsub run-visits.sh -u 2.56e-01 --effects=1; done
for i in `seq 1 10`; do qsub run-visits.sh -u 5.0e-02 --effects=1; done
