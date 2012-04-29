# this estimates the probability distribution of a sample x where x is a non-negative discrete random variable
comp.distr <- function(x) {
  x.distr <- numeric(max(x)+1)
  tbl <- table(x)
  x.distr[as.numeric(names(tbl))+1] <- as.numeric(tbl)
  x.distr/sum(x.distr)
}

library(grid)

N <- 1000

pip <- pipe("ls out/sojourn", "r")
sets <- scan(file=pip, what="", sep="\n")
close(pip)

mutation.rates <- sub(".*-u-", "", sets)

sets.data <- lapply(sets, function(set) {
  cat("processing set", set, "\n")
  set.dir <- paste("out/sojourn/", set, sep="")
  
  # get file listings
  pip <- pipe(paste("lfs find", set.dir, "-name 'equilibrium*.out.gz'"), "r")
  files <- scan(file=pip, what="", sep="\n")
  close(pip)
  sojourn.files <- sub("out.gz$", "sojourns.gz", files)
  segsite.files <- sub("out.gz$", "segsites.gz", files)
  
  # read in the sojourn tables
  sojourn.d <- lapply(sojourn.files, function(f) read.table(file=f, header=FALSE))
  sojourn.d <- do.call(rbind, sojourn.d)
  names(sojourn.d) <- c("time", "effect")
  max.sojourn <- max(sojourn.d$time)
  num.sojourns <- nrow(sojourn.d)
  
  # Currently I assume there is one effect size
  stopifnot(length(unique(abs(sojourn.d$effect))) == 1)

  # create a table of sojourn distributions
  sojourn.distr <- comp.distr(sojourn.d$time)

  # read in segregating sites tables
  segsites <- unlist(lapply(segsite.files, function(f) scan(file=f, what=0, sep="\n")))
  observations <- length(segsites)
  segsites.distr <- comp.distr(segsites)
  
  list(sojourn.distr=sojourn.distr, segsites.distr=segsites.distr, observations=observations, num.sojourns=num.sojourns)
})

sampled.segsite.distr <- lapply(1:length(sets.data), function(i) {
  cat("sampling for dataset", i, "\n")
  N <- 1000
  d <- sets.data[[i]]
  mu <- as.numeric(mutation.rates)[i]
  rcdf <- rev(cumsum(rev(d$sojourn.distr)))
  s <- sapply(1:d$observations, function(i) sum(rpois(length(rcdf)-1, lambda=2*N*mu*rcdf[-1])))
  comp.distr(s)
})

neg.inf.val <- -7
palette(c("purple", "tan", "olivedrab", "gray30"))
pdf(file="seg_sites_resampling_from_sojourn_times.pdf", height=4, width=4.5)
par(mar=c(4,4,1,4), mgp=c(2,0.8,0), cex.lab=0.95, cex.axis=0.75, xpd=NA)
plot(c(neg.inf.val-1,0), c(neg.inf.val-1,0), type="n", axes=FALSE,
  xlab="log10 empirical distribution", ylab="log10 resampled distribution")
wd <- 0.3
# vertical gray bar highlighting the -Inf region
rect(neg.inf.val-wd, neg.inf.val-2.5, neg.inf.val+wd, 0, border="gray90", col="gray90")
# horizontal gray bar highlighting the -Inf region
rect(neg.inf.val-2.1, neg.inf.val-wd, 0, neg.inf.val+wd, border="gray90", col="gray90")
segments(neg.inf.val, neg.inf.val, 0, 0, lwd=1, lty=2)
axis.at <- seq(neg.inf.val-1, 0, by=2)
axis.at[1] <- neg.inf.val
axis.labs <- axis.at
axis.labs[1] <- expression(-infinity)
axis(side=1, at=axis.at, labels=axis.labs) 
axis(side=2, at=axis.at, labels=axis.labs) 
for (i in length(sets.data):1) {
  # select the data we want
  x1 <- sets.data[[i]]$segsites.distr
  y1 <- sampled.segsite.distr[[i]]

  # make the distribution vectors the same link (we know they both start at zero)
  len <- max(length(x1),length(y1))
  x2 <- numeric(len)
  y2 <- numeric(len)
  x2[1:length(x1)] <- x1
  y2[1:length(y1)] <- y1

  # take the log10 and set the -Inf to a special value
  x2 <- log10(x2)
  y2 <- log10(y2)
  stopifnot(min(x2[x2>-Inf]) > neg.inf.val)
  stopifnot(min(y2[y2>-Inf]) > neg.inf.val)
  # I don't plot points where both are -Inf, but cases where only one is -Inf, I use the special value
  sel <- x2 == -Inf & y2 == -Inf
  x2[x2 == -Inf] <- neg.inf.val
  y2[y2 == -Inf] <- neg.inf.val
  x2[sel] <- NA
  y2[sel] <- NA

  # plot it
  lines(x2, y2, col=i, lwd=2)

  # plot a black dot at the beginning of the distribution (the lowest number of segregating sites
  points(x2[!sel][1], y2[!sel][1], pch=20, col=i, cex=1.1)
}
o <- order(as.numeric(mutation.rates))
legend(bty="n", cex=0.7, legend=paste("2Nu =", 2*N*as.numeric(mutation.rates)[o]), 
  col=(1:length(sets.data))[o], lwd=2, x=2, y=neg.inf.val+0.5, xjust=1, yjust=0)
dev.off()

# prepare to plot the log densities individually
max.segsites.empirical <- sapply(sets.data, function(x) length(x$segsites.distr))-1
max.segsites.resampled <- sapply(sampled.segsite.distr, length)-1
max.segsites <- max(c(max.segsites.empirical, max.segsites.resampled))

# find minimum non-infinite densities 
min.finitie.log.density.empirical <- min(sapply(sets.data, function(x) {
  dens <- log10(x$segsites.distr)
  min(dens[dens>-Inf])
}))
min.finitie.log.density.resampled <- min(sapply(sampled.segsite.distr, function(x) {
  dens <- log10(x)
  min(dens[dens>-Inf])
}))
min.finitie.log.density <- min(min.finitie.log.density.empirical, min.finitie.log.density.resampled)

# plot the log densities individually
pdf(file="seg_sites_resampling_from_sojourn_times-individual.pdf", height=3, width=5)
par(mar=c(4,4,1,4), mgp=c(2,0.8,0), cex.lab=0.95, cex.axis=0.75, xpd=NA)
plot(c(0, max.segsites), c(min.finitie.log.density,0), type="n", bty="n",
  xlab="number of segregating sites", ylab="log10 density")
for (i in 1:length(sets.data)) {
  dens <- sets.data[[i]]$segsites.distr
  lines(0:(length(dens)-1), log10(dens), col=i, lty=1)
  dens <- sampled.segsite.distr[[i]]
  lines(0:(length(dens)-1), log10(dens), col=i, lty=2)
}
o <- order(as.numeric(mutation.rates))
legend(x=750, y=0.1, xjust=0, yjust=1, bty="n", cex=0.6, legend=paste("2Nu =", 2*N*as.numeric(mutation.rates)[o]), 
  col=(1:length(sets.data))[o], lwd=1)
legend(x=790, y=0.1, xjust=1, yjust=1, bty="n", cex=0.6, legend=rep("", length(mutation.rates)), 
  col=(1:length(sets.data))[o], lwd=1, lty=2)
text("resampled", x=700, y=-0.1, cex=0.6, srt=30, adj=c(0,0))
text("empirical", x=810, y=-0.1, cex=0.6, srt=30, adj=c(0,0))
dev.off()

# END


