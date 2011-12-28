library(grid)

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
  
  # read in the sojourn tables
  sojourn.d <- lapply(sojourn.files, function(f) read.table(file=f, header=FALSE))
  sojourn.d <- do.call(rbind, sojourn.d)
  names(sojourn.d) <- c("time", "effect")
  
  # use the absolute value of the effect size to split up sojourn times
  sojourn.by.effect <- split(sojourn.d$time, abs(sojourn.d$effect))
  
  max.sojourn <- max(sojourn.d$time)
  
  # create a table of sojourn distributions
  sojourn.distr <- matrix(0, nrow=max.sojourn, ncol=length(sojourn.by.effect))
  for (i in 1:length(sojourn.by.effect)) {
    tbl <- table(sojourn.by.effect[[i]])
    sojourn.distr[as.numeric(names(tbl)),i] <- as.numeric(tbl)
  }
  sojourn.distr <- apply(sojourn.distr, 2, function(x) x/sum(x))
  
  sojourn.distr.log <- log10(sojourn.distr)
  min.nonzero.density.log <- min(sojourn.distr.log[sojourn.distr.log > -Inf])
  max.density.log <- max(sojourn.distr.log)
  list(distr=sojourn.distr, min.density=min.nonzero.density.log, max.density=max.density.log, max.sojourn=max.sojourn)
})

max.density.log <- max(sapply(sets.data, function(x) x$max.density))
min.nonzero.density.log <- min(sapply(sets.data, function(x) x$min.density))
max.sojourn <- max(sapply(sets.data, function(x) x$max.sojourn))

# set up a pretty log-scale axis
y.scale <- c(min.nonzero.density.log,max.density.log)
axis.ticks <- t(matrix(10^pretty(y.scale)[-1], ncol=1) %*% matrix(seq(0,1,length.out=11)[-1], nrow=1))
axis.ticks <- sort(unique(axis.ticks[1:prod(dim(axis.ticks))]))
axis.ticks.lab <- format(axis.ticks)
axis.ticks.lab[round(log10(axis.ticks)) != log10(axis.ticks)] <- NA

# plot the sojourn distributions
palette(c("gray30", "purple", "firebrick", "dodgerblue"))
scale <- 1.3
pdf(file="sojourn_distributions.pdf", height=4*scale, width=7*scale)
pushViewport(viewport(x=0.97, y=0.97, just=c(1,1), height=0.8, width=0.8, 
    xscale=c(0,max.sojourn), yscale=range(pretty(y.scale))))
  trash <- lapply(1:length(sets.data), function(i) {
    grid.points(x=1:nrow(sets.data[[i]]$distr), y=log10(sets.data[[i]]$distr[,1]), 
      default.units="native", size=unit(0.3, "char"), pch=21, gp=gpar(col=i))
  })
#  grid.text(paste("mutations:", format(as.real(nrow(sojourn.d)), scientific=TRUE, digits=4)), 
#    x=1, y=1, just=c(1,1), gp=gpar(cex=0.8))
  grid.xaxis(gp=gpar(cex=0.8))
  grid.yaxis(at=log10(axis.ticks), label=FALSE, gp=gpar(cex=0.8))
  grid.text(axis.ticks.lab[axis.ticks.lab!=""], y=log10(axis.ticks)[axis.ticks.lab!=""], 
    x=unit(-0.06, "npc"), gp=gpar(cex=0.8), default.units="native")
  grid.text("sojourn time (generations)", x=0.5, y=-0.15, gp=gpar(cex=1.0))
  grid.text("density", y=0.5, x=-0.15, rot=90, gp=gpar(cex=1.0))
popViewport()
# plot a legend
pushViewport(viewport(x=0.97, y=0.97, just=c(1,1), height=0.2, width=0.2))
  mutation.o <- order(as.numeric(mutation.rates))
  legend.y <- seq(0, 0.7, length.out=length(sets))
  grid.points(x=rep(0, length(sets)), y=legend.y, size=unit(0.4, "char"), pch=21, gp=gpar(col=(1:length(sets))[mutation.o]))
  grid.text(mutation.rates[mutation.o], x=0.1, just=c(0,0.5), y=legend.y, gp=gpar(cex=0.8, fontfamily="mono"))
  grid.text("mutation rate", x=0, y=1, just=c(0,0.5))
popViewport()
dev.off()

# END


