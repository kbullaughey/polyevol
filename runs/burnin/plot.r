library(grid)
source("../../src/r/evolveq.R")

# read in output files
pip <- pipe("find . -name 'burnin*.out'", "r")
files <- scan(pip, what="", sep="\n")
close(pip)

runs <- lapply(files, load.simoutput)
nloci <- runs[[1]]$par$loci
sum.het <- numeric(runs[[1]]$par$times+1)

for (i in 1:length(runs)) {
  for (l in 1:nloci) {
    p <- runs[[i]]$freq[[l]]$freq
    # for the generations that this site is segregating, add the heterozygosity to sum.het 
    sum.het[runs[[i]]$freq[[l]]$gen+1] <- sum.het[runs[[i]]$freq[[l]]$gen+1] + 2*p*(1-p)
  }
}

# END
