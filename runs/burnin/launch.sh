for i in `seq 1 10`; do for s in `cat s_values`; do qsub burnin_run.sh $s 10; done; done
