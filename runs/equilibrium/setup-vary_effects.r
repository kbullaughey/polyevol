# Return a modified model, setting mu so that the expected mutational input 
# averaging over effects is the mu given in the model, m.
adjust.mu <- function(m) {
  mu <- signif(m$mu*1/sum(m$effects * m$eprobs), digits=4)
  list(effects=m$effects, eprobs=m$eprobs, mu=mu)
}

# Print the parameters included in m, in the command line format
format.params <- function(m) {
  p <- sapply(m, paste, collapse=",", sep="")
  paste("--", names(m), "=", p, sep="", collapse=" ")
}

# combinations of effect sizes
effects <- list(
  1,
  2,
  0.5,
  0.1,
  c(1,0.1),
  c(1,0.2),
  c(2,0.1),
  c(2,0.2)
)

# combinations of effect-size probability distributions
eprobs <- list(
  1,
  c(0.9,0.1),
  c(0.5,0.5),
  c(0.1,0.9)
)

# vector of mutation rates to choose among
mu.vec <- c(0.1, 0.5, 1, 2, 10)/2000

# generate models for all combinations of effects, eprobs, and mu
g <- expand.grid(1:length(effects), 1:length(eprobs), 1:length(mu.vec))
param.sets <- lapply(1:nrow(g), function(i) {
  list(effects=effects[[g[i,1]]], eprobs=eprobs[[g[i,2]]], mu=mu.vec[g[i,3]])
})

# limit the models only to those that have compatable effects and eprobs vectors
effects.len <- sapply(param.sets, function(x) length(x$effects))
eprobs.len <- sapply(param.sets, function(x) length(x$eprobs))
param.sets <- param.sets[effects.len == eprobs.len]

# adjust mu in each param set so that the expected mutational input (averaging effects) is the given mu
param.sets.adj <- lapply(param.sets, adjust.mu)
total.mu <- round(sapply(param.sets.adj, function(x) sum(x$effects * x$eprobs * x$mu)), 2)

# write a file with all the commands to lanuch these jobs
cat(file="launch-vary_effects.sh", paste("qsub run-vary_effects.sh", sapply(param.sets.adj, format.params), "\n"), sep="")


# END
