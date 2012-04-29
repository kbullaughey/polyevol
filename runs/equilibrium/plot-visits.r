library(grid)
source('../../src/r/evolveq.R')

run.name <- 'visits-2012_01_31'
pip <- pipe(paste("ls out/", run.name, sep=""), "r")
sets <- scan(file=pip, what="", sep="\n")
close(pip)

mutation.rates <- as.numeric(sub(".*u-", "", sapply(strsplit(sets, split="--"), 
  function(x) grep("u-", x, value=TRUE))))

sets.data <- lapply(sets, function(set) {
  cat("processing set", set, "\n")
  set.dir <- paste("out/", run.name, "/", set, sep="")
  
  # get file listings
  # TODO: change find to lfs find
  pip <- pipe(paste("find", set.dir, "-name 'visits*.out.gz'"), "r")
  files <- scan(file=pip, what="", sep="\n")
  close(pip)

  # read in all the files for this set
  runs.lines <- lapply(files, function(f) scan(file=f, what="", sep="\n"))
  runs.params <- lapply(runs.lines, parse.params)

  # get the mutation rate and make sure all runs are the same
  mu <- unique(sapply(runs.params, function(x) x$mu))
  stopifnot(length(mu) == 1)

  # make sure population sizes are the same
  popsize <- unique(sapply(runs.params, function(x) x$popsize))
  
  # extract visit vectors
  visit.lines <- grep("visits:", unlist(runs.lines), value=TRUE)
  visit.data <- lapply(strsplit(sub(".*visits: ", "", visit.lines), split=" "), as.numeric)
  stopifnot(sum(sapply(visit.data, length) != 2*popsize-1) == 0)
  visit.data <- t(matrix(unlist(visit.data), ncol=2*popsize-1, byrow=TRUE))

  # compute overall visits distribition
  visits <- apply(visit.data, 1, sum)
  visits.distr <- visits / (mu * 2*popsize * runs.params[[1]]$times[1] * length(files))
  visits.distr.log <- log10(visits.distr)
  min.nonzero.density.log <- min(visits.distr.log[visits.distr.log > -Inf])
  max.density.log <- max(visits.distr.log)
  
  # read in the sojourn tables
  list(distr=visits.distr.log, min.density=min.nonzero.density.log, max.density=max.density.log, 
    mu=mu, N=popsize, params=runs.params)
})

N <- sapply(sets.data, function(x) x$N)
mu <- sapply(sets.data, function(x) x$mu)
max.density.log <- max(sapply(sets.data, function(x) x$max.density))
min.nonzero.density.log <- min(sapply(sets.data, function(x) x$min.density))

# so far I only support one population size
N <- unique(N)
stopifnot(length(N) == 1)

# set up a pretty log-scale axis
y.scale <- c(min.nonzero.density.log,max.density.log)
axis.ticks <- t(matrix(10^seq(min(pretty(y.scale)), max(pretty(y.scale)))[-1], ncol=1) %*% 
  matrix(seq(0,1,length.out=11)[-1], nrow=1))
axis.ticks <- sort(unique(signif(axis.ticks[1:prod(dim(axis.ticks))], digits=1)))
axis.ticks.lab <- format(axis.ticks)
axis.ticks.lab[floor(log10(axis.ticks)) != log10(axis.ticks)] <- NA

# plot the sojourn distributions
palette(rep(c("gray30", "purple", "firebrick", "dodgerblue", "tan", "orange", "darkblue", "violet", "olivedrab"), 10))
scale <- 1.3
pdf(file=paste("visits_distributions-", run.name, ".pdf", sep=""), height=4*scale, width=7*scale)
pushViewport(viewport(x=0.97, y=0.97, just=c(1,1), height=0.8, width=0.8, 
    xscale=c(0,2*N), yscale=range(pretty(y.scale))))
  derived.allele.counts <- 1:(2*N-1)
  trash <- lapply(1:length(sets.data), function(i) {
    grid.lines(x=derived.allele.counts, y=sets.data[[i]]$distr, 
      default.units="native", gp=gpar(col=i, lwd=1.5))
  })
  neutral.expectation <- 2/derived.allele.counts
  neutral.expectation <- neutral.expectation / sum(neutral.expectation) * log(2*N)
  grid.lines(x=1:(2*N-1), y=log10(neutral.expectation), default.units="native", gp=gpar(col="black", lwd=2, lty=2))
  grid.xaxis(gp=gpar(cex=0.8))
  grid.yaxis(at=log10(axis.ticks), label=FALSE, gp=gpar(cex=0.8))
  grid.text(axis.ticks.lab[axis.ticks.lab!=""], y=log10(axis.ticks)[axis.ticks.lab!=""], 
    x=unit(-0.06, "npc"), gp=gpar(cex=0.8), default.units="native")
  grid.text("Derived allele count (freq. bin)", x=0.5, y=-0.15, gp=gpar(cex=1.0))
  grid.text("density", y=0.5, x=-0.15, rot=90, gp=gpar(cex=1.0))
popViewport()
# plot a legend
pushViewport(viewport(x=0.90, y=0.99, just=c(1,1), height=0.45, width=0.1))
  mutation.o <- order(as.numeric(mutation.rates * 2*N))
  legend.y <- seq(0, 0.75, length.out=length(sets))
  grid.points(x=rep(0.2, length(sets)), y=legend.y, size=unit(0.4, "char"), pch=21, gp=gpar(col=(1:length(sets))[mutation.o]))
  grid.text(2*N*mutation.rates[mutation.o], x=0.3, just=c(0,0.5), y=legend.y, gp=gpar(cex=0.8, fontfamily="mono"))
  grid.text("mutations entering", x=0.5, y=1, just=c(0.5,1), gp=gpar(cex=0.9))
  grid.text("population per gen.", x=0.5, y=0.9, just=c(0.5,1), gp=gpar(cex=0.9))
popViewport()
dev.off()

# END


