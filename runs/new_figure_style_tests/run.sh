# Run four sets of 5 runs
# In order for a polymorphism not to end up balanced, the optimum must change 
#   by 2, when it changes by one there will be a polymorphsm balanced before 
#   or after the change in the optimum (depending on whether the optimum is 
#   even or odd and whether it's the infinite sites or finite sites model.

# infinite sites, a change to an optimum that will result in a balanced polymorphism
#for i in `seq 1 5`; do 
#  echo $i
#  quant --seed=`rand.pl` -m infinite -N 4000 -s 0.07 --effects=1 --loci=25 \
#    --freqs=even --opts=0,1 --times=300,1500 --burnin=2000 \
#   > balancing-infinite.$i.out
#done
#
## finite sites, a change to an optimum that will result in a balanced polymorphism
#for i in `seq 1 5`; do 
#  echo $i
#  quant --seed=`rand.pl` -m finite -N 4000 -s 0.07 --effects=1 --loci=25 \
#    --freqs=even --opts=1,2 --times=300,1500 --burnin=2000 \
#   > balancing-finite.$i.out
#done

# infinite sites, no balanced polymorphisms, a larger change in the optimum
for i in `seq 1 5`; do 
  echo $i
  quant --seed=`rand.pl` -m infinite -N 4000 -s 0.07 --effects=1 --loci=25 \
    --freqs=even --opts=0,2 --times=300,1500 --burnin=2000 \
   > no_balancing-infinite.$i.out
done

# finite sites, no balanced polymorphisms, a larger change in the optimum
for i in `seq 1 5`; do 
  echo $i
  quant --seed=`rand.pl` -m finite -N 4000 -s 0.07 --effects=1 --loci=25 \
    --freqs=even --opts=1,3 --times=300,1500 --burnin=2000 \
   > no_balancing-finite.$i.out
done


# END
