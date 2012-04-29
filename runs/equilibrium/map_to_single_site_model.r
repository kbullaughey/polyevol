library(grid)
source('../../src/r/evolveq.R')

run.name <- 'combined_visits'
pip <- pipe("(lfs find out/visits; lfs find out/visits-2012_01_31) | grep effects_1 | grep 'out.gz$'", "r")
files <- scan(file=pip, what="", sep="\n")
close(pip)

mutation.rates <- as.numeric(sub("^.*visits-[-_0-9]*u-(.*)--effects_1\\.[0-9]*\\.out\\.gz$", "\\1", files))
files.split <- split(files, mutation.rates)
set.mu <- as.numeric(names(files.split))

sets.data <- lapply(files.split, function(files) {
  # read in all the files for this set
  runs.lines <- lapply(files, function(f) scan(file=f, what="", sep="\n"))
  runs.params <- lapply(runs.lines, parse.params)

  # get the mutation rate and make sure all runs are the same
  mu <- unique(sapply(runs.params, function(x) x$mu))
  stopifnot(length(mu) == 1)

  # make sure population sizes are the same
  popsize <- unique(sapply(runs.params, function(x) x$popsize))
  stopifnot(length(popsize) == 1)
  
  # extract visit vectors
  visit.lines <- grep("visits:", unlist(runs.lines), value=TRUE)
  visit.data <- lapply(strsplit(sub(".*visits: ", "", visit.lines), split=" "), as.numeric)
  stopifnot(sum(sapply(visit.data, length) != 2*popsize-1) == 0)
  visit.data <- t(matrix(unlist(visit.data), ncol=2*popsize-1, byrow=TRUE))

  # compute overall visits distribition
  visits <- apply(visit.data, 1, sum)
  
  # read in the sojourn tables
  list(visits=visits, mu=mu, N=popsize, params=runs.params)
})

N <- sapply(sets.data, function(x) x$N)
mu <- sapply(sets.data, function(x) x$mu)

# make sure mutation rates all match
stopifnot(sum(mu != set.mu) == 0)

# so far I only support one population size
N <- unique(N)
stopifnot(length(N) == 1)

site.density <- function(gamma, q) 
  (1 - exp(-2*gamma*(1-q))) / ((1 - exp(-2*gamma)) * q * (1-q))

# Approxmiation given by Srinivasa Ramanujan. Source: http://en.wikipedia.org/wiki/Factorial
log.fac.approx <- function(n)
  n * log(n) - n + log(n*(1+4*n*(1+2*n)))/6 + log(pi)/2

# For comparison. 
log.fac.exact <- function(n) {
  sum(sapply(1:n, log))
}

# Prob(visits | gamma)
log.likelihood <- function(gamma, visits) {
  i <- 1:(2*N-1)
  n.i <- visits
  stopifnot(length(i) == length(visits))
  p.i <- sapply(i/(2*N), function(q) site.density(gamma, q))
  p.i <- p.i / sum(p.i)
  nonzero <- visits > 0
  sum(n.i[nonzero] * log(p.i[nonzero]) - log.fac.approx(n.i[nonzero]) - p.i[nonzero])
}

# compute likelihood surfaces
gamma.grid <- seq(-200, -1, by=1)
surfaces <- lapply(sets.data, function(x) {
  cat(".")
  sapply(gamma.grid, log.likelihood, x$visits)
})
cat("\n")

# compute MLE estimates for the PRF gamma that corresponds to each polygenic model
mle <- sapply(sets.data, function(x) 
  optimize(log.likelihood, interval=c(-200, -0.01), x$visits, maximum=TRUE)$maximum)

save.image(file="map_to_single_site_model.rimage")

# look at the likelihood surface
palette(topo.colors(length(surfaces)))
pdf(file="likelihood_surfaces.pdf", height=3, width=4)
par(mar=c(4,4,1,1), mgp=c(2,0.8,0), cex.lab=1.1, cex.axis=0.8)
trash <- lapply(1:length(surfaces), function(i) {
  plot(range(gamma.grid), range(surfaces[[i]])/10^6, type="n", 
    xlab="gamma", ylab="log likelihood (10^6)")
  lines(gamma.grid, surfaces[[i]]/10^6, col="black")
})
dev.off()

# plot mutational input vs gamma
pdf(file="mutational_input_vs_gamma.pdf", height=4, width=4)
par(mar=c(4,4,1,1), mgp=c(2,0.8,0), cex.lab=1.0, cex.axis=0.75)
plot(set.mu * 2*N, mle, pch=21, cex=0.8, col="purple", 
  xlab=expression(paste(2, plain(N), mu)), ylab=expression(paste(2, N, gamma)))
dev.off()

# look at the log site density
gamma.grid <- c(-0.08, -0.04, -0.02, -0.01, -0.005, -0.0025, -0.00125) * 2 * N
freq.bins <- 1:(2*N-1)
log.site.densities <- lapply(gamma.grid, function(gam) {
  log10(site.density(gam, freq.bins/(2*N)))
})
pdf(file="log_site_density.pdf", height=3, width=5)
par(mar=c(4,4,1,1), mgp=c(2,0.8,0), cex.lab=1.1, cex.axis=0.8)
plot(c(0,2*N), c(-12,3), type="n", lwd=2, col="purple", 
  xlab="derived allele count", ylab="log density")
trash <- lapply(log.site.densities, function(log.site.density) {
  sel <- log.site.density > -12
  lines(freq.bins[sel], log.site.density[sel], lwd=1, col="black")
})
dev.off()

# look at the first few log visits densities
pdf(file="log_visits_densities.pdf", height=3, width=5)
par(mar=c(4,4,1,1), mgp=c(2,0.8,0), cex.lab=1.1, cex.axis=0.8)
plot(c(0,2*N), c(-12,1), type="n", lwd=2, col="purple", 
  xlab="derived allele count", ylab="log density")
trash <- lapply(sets.data[1:10], function(x) {
  dens <- log10(x$visits / sum(x$visits))
  sel <- dens > -12
  lines(freq.bins[sel], dens[sel], lwd=1, col="black")
})
dev.off()

# Compare empirical distribution and several theoretical distributions
neut <- 2/freq.bins
neut <- neut/sum(neut)
g.n20 <- site.density(-20, freq.bins/(2*N))
g.n20 <- g.n20/sum(g.n20)
g.n100 <- site.density(-100, freq.bins/(2*N))
g.n100 <- g.n100/sum(g.n100)
emp <- sets.data[[1]]$visits
emp <- emp/sum(emp)

palette(c("purple", "tan", "olivedrab"))
pdf(file="distr_compare.pdf", height=4, width=4)
par(mar=c(4,4,1,1), mgp=c(2,0.8,0), cex.lab=1.0, cex.axis=0.75)
plot(c(-10,0), c(-10,0), type="n", xlab="log10 empirical distribution", ylab="log10 PRF distribution")
abline(0,1,col="black", lwd=1, lty=2)
lines(log10(emp), log10(neut), col=1, lwd=2)
lines(log10(emp), log10(g.n20), col=2, lwd=2)
lines(log10(emp), log10(g.n100), col=3, lwd=2)
legend("bottomright", inset=0.03, bty="n", 
  legend=c(expression(plain(neutral)), expression(gamma == -20), expression(gamma == -100)), 
  col=1:3, lwd=2)
dev.off()

emp.all <- lapply(sets.data, function(z) {
  z$visits / sum(z$visits)
})
max.emp.pretty <- max(pretty(sapply(emp.all, max), 2))


# Compare each empirical distribution to the best fit
palette(c("olivedrab", rgb(0.97, 0.97, 0.97)))
scale <- 0.8
pdf(file="distr_compare-all.pdf", height=3.1*scale, width=9*scale)
par(mar=c(4,4,2.2,1), mgp=c(2,0.8,0), cex.lab=1.0, cex.axis=0.75, mfcol=c(1,3))
for (i in 1:length(sets.data)) {
  gamma.mle.distr <- site.density(mle[[i]], freq.bins/(2*N))
  gamma.mle.distr <- gamma.mle.distr / sum(gamma.mle.distr)
  emp <- emp.all[[i]]

  # plot log scale
  high <- log10(max.emp.pretty)
  plot(c(-10,0), c(-10,0), type="n", xlab="log10 empirical distribution", ylab="log10 PRF distribution")
  rect(-3, -3, high, high, col=2, border=NA)
  segments(c(-3, -3, -3, high), c(-3, -3, high, high), c(-3, high, high, high), c(high, -3, high, -3), 
    lty=2, col="gray")
  abline(0,1,col="black", lwd=1, lty=2)
  lines(log10(emp), log10(gamma.mle.distr), col=1, lwd=2)
  legend("bottomright", inset=0.03, bty="n", 
    legend=substitute(gamma == m, list(m=signif(mle[[i]], digits=3))), 
    col=1, lwd=2)

  # plot regular scale
  high <- max.emp.pretty
  plot(c(0.001,high), c(0.001,high), type="n", xlab="empirical distribution", ylab="PRF distribution", 
    main=paste("2Nu =", as.numeric(set.mu[i])*2*N), axes=FALSE)
  rect(-0.1, -0.1, 1, 1, col=2, border=NA)
  box()
  axis(1)
  axis(2)
  abline(0,1,col="black", lwd=1, lty=2)
  sel <- 0.001 < emp & emp < high & 0.001 < gamma.mle.distr & gamma.mle.distr < high
  lines(emp[sel], gamma.mle.distr[sel], col=1, lwd=2)
  legend("bottomright", inset=0.03, bty="n", 
    legend=substitute(gamma == m, list(m=signif(mle[[i]], digits=3))), 
    col=1, lwd=2)

  # plot the likelihood surface
  plot(range(gamma.grid), range(surfaces[[i]])/10^6, type="n", 
    xlab="gamma", ylab="log likelihood (x10^6)")
  lines(gamma.grid, surfaces[[i]]/10^6, col="black")
  abline(v=mle[[i]], col="gray", lty=2)
}
dev.off()

# END


