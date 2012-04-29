source("../../src/r/evolveq.R")

# combinations of effect sizes
effects <- seq(0, 2, by=0.05)[-1]

# combinations of mutation rates
mu <- seq(0, 0.16, by=0.004)[-1]


# generate models for all combinations of effects, eprobs, and mu
g <- expand.grid(1:length(effects), 1:length(mu))
param.sets <- lapply(1:nrow(g), function(i) {
  list(effects=effects[[g[i,1]]], mu=mu[g[i,2]])
})

# write a file with all the commands to lanuch these jobs
cat(file="launch-phenotype_effect_grid.sh", paste("qsub run-phenotype_effect_grid.sh", 
  sapply(param.sets, format.params), "\n"), sep="")

save.image(file="setup-phenotype_effect_grid.rimage")


# END
