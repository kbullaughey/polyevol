# This script produces a plot of how the number of segregating sites scales with the mutation rate for several effect sizes including 0, which is equivalent to neutrality


# load necessary libraries
source("../../src/r/evolveq.R")
library(grid)

# read in the file names for all simulation output files
pip <- pipe("lfs find out/segregating_sites -name 'segsite_counts-*.out.gz'", "r")
files <- scan(file=pip, what="", sep="\n")
close(pip)

# read in each file
data <- lapply(files, function(f) {
  file.lines <- readLines(f)
  params <- parse.params(file.lines)
  segregating.sites <- as.numeric(file.lines[-(1:2)])
  list(params=params, mean=mean(segregating.sites), sd=sd(segregating.sites))
})

# split the list by effect size and then sort by mutation rate, so points 
#   will be in the correct order in the plot
data.split <- split(data, sapply(data, function(x) x$params$effects))
data <- lapply(data.split, function(d) {
  d[order(sapply(d, function(x) x$params$mu))]
})

# extract the data we care about into matricies, with one column per effect size
N <- sapply(data, function(d) sapply(d, function(x) x$params$popsize))
u <- sapply(data, function(d) sapply(d, function(x) x$params$mu))
ss.mean <- sapply(data, function(d) sapply(d, function(x) x$mean))
ss.sd <- sapply(data, function(d) sapply(d, function(x) x$sd))

# for sapply to build matricies, each effect size must have the same length
stopifnot(length(dim(N))==2)

# check that we have the same selection of simulations for each effect size
stopifnot(sum(sapply(2:ncol(N), function(i) sum(N[,i]!=N[,1]))) == 0)
stopifnot(sum(sapply(2:ncol(u), function(i) sum(u[,i]!=u[,1]))) == 0)

theta <- 4*u[,1]*N[,1]
neutral.ev <- theta * (log(2*N[,1]) + 0.6775)

# plot the number of segregating sites as a function of theta
palette(c("black", "gray90", "gray60", "dodgerblue", "orange", "tan4"))
pdf(file="segregating_sites-mutation_rate.pdf", width=6, height=4)
pushViewport(viewport(x=0.98, y=0.98, just=c(1,1), height=0.8, width=0.8, 
    xscale=c(0,max(theta)), yscale=c(0,max(ss.mean+ss.sd))))
  # plot the neutral expectation
  grid.polygon(c(theta, rev(theta)), c(neutral.ev+sqrt(neutral.ev), rev(neutral.ev-sqrt(neutral.ev))), 
    default.units="native", gp=gpar(col=4, fill=3))
  grid.lines(theta, neutral.ev, default.units="native", gp=gpar(col=2))

  # plot the simulation results
  trash <- lapply(1:ncol(N), function(i) {
    grid.segments(theta, ss.mean[,i]+ss.sd[,i], theta, ss.mean[,i]-ss.sd[,i], 
      gp=gpar(col=1), default.units="native")
    grid.points(theta, ss.mean[,i], default.units="native", pch=20, size=unit(0.6, "char"), 
      gp=gpar(col=1))
  })
  grid.xaxis(gp=gpar(cex=0.8))
  grid.yaxis(gp=gpar(cex=0.8))
  grid.text("segregating sites", y=0.5, x=-1.1, rot=90)
  grid.text(expression(theta = Nu), x=0.5, y=-1.1)
popViewport()
dev.off()


# END
