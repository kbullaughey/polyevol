library(grid)
source("../../src/r/evolveq.R")

the.args <- commandArgs()
model <- sub("--model=", "", the.args[grep("--model=", the.args)])

# read in output files
pip <- pipe(paste("find out-", model, " -name 'burnin*.out'", sep=""), "r")
files <- scan(pip, what="", sep="\n")
close(pip)

runs <- lapply(files, load.simoutput)
save.image(file=paste("burnin_plot_data-", model, ".rimage", sep=""))

nloci <- runs[[1]]$par$loci
generations <- runs[[1]]$par$times
file.s <- as.numeric(sub("^.*out-[a-z]*/([.0-9]*)/.*$", "\\1", files))
runs.split <- split(runs, file.s)
num.s <- length(runs.split)

sum.hets <- lapply(runs.split, function(run.set) {
  sum.het <- numeric(generations+1)
  for (i in 1:length(run.set)) {
    for (l in 1:nloci) {
      p <- run.set[[i]]$freq[[l]]$freq
      # for the generations that this site is segregating, add the heterozygosity to sum.het 
      sum.het[run.set[[i]]$freq[[l]]$gen+1] <- sum.het[run.set[[i]]$freq[[l]]$gen+1] + 2*p*(1-p)
    }
  }
  return(sum.het/length(run.set))
})

max.het <- max(unlist(sum.hets))

# get some colors, throwing away the first and last 
palette(terrain.colors(num.s+2)[1:num.s])

# plot sum (over loci) heterozygosity as a function of generation for each s
pdf(file=paste("burnin_sum_heterozygosity-", model, ".pdf", sep=""), width=6, height=4)
pushViewport(viewport())
  pushViewport(viewport(x=0.95, y=0.95, just=c(1,1), height=0.75, width=0.8))
    pushViewport(viewport(xscale=c(0,generations), yscale=c(0,max.het)))
      for (i in 1:num.s) {
        grid.lines(x=0:generations, y=sum.hets[[i]], gp=gpar(col=i, lwd=2), default.units="native")
      }
      grid.xaxis(gp=gpar(cex=0.8))
      grid.yaxis(gp=gpar(cex=0.8))
    popViewport()
    grid.text("generations", x=0.5, y=-0.18)
    grid.text("avg. sum (over loci) heterozygosity", y=0.5, x=-0.13, rot=90)
  popViewport()
  grid.text(expression(paste("selection: ", Ns, alpha^2)), x=0.825, y=0.92, gp=gpar(cex=0.9))
  pushViewport(viewport(x=0.95, y=0.9, just=c(1,1), height=0.2, width=0.2))
    legend.y <- (1:num.s)/(num.s+1)
    grid.segments(rep(0.1, num.s), legend.y, rep(0.3, num.s), legend.y, gp=gpar(col=1:num.s, lwd=2))
    s <- as.numeric(names(runs.split))
    Nsa2 <- runs[[1]]$par$popsize * s * runs[[1]]$par$effects^2
    grid.text(Nsa2, x=0.4, y=legend.y, just=c(0,0.5), gp=gpar(cex=0.75))
  popViewport()
popViewport()
dev.off()

# plot some example trajecories (not averaged)

sum.hets.indiv <- lapply(runs.split, function(run.set) {
  lapply(run.set, function(run) {
    sum.het <- numeric(generations+1)
    for (l in 1:nloci) {
      p <- run$freq[[l]]$freq
      # for the generations that this site is segregating, add the heterozygosity to sum.het 
      sum.het[run$freq[[l]]$gen+1] <- sum.het[run$freq[[l]]$gen+1] + 2*p*(1-p)
    }
    return(sum.het)
  })
})
max.het.indiv <- max(unlist(sum.hets.indiv))

pdf(file=paste("burnin_sum_heterozygosity-individual-", model, ".pdf", sep=""), width=6, height=4)
pushViewport(viewport())
  pushViewport(viewport(x=0.95, y=0.95, just=c(1,1), height=0.75, width=0.8))
    pushViewport(viewport(xscale=c(0,generations), yscale=c(0,max.het.indiv)))
      for (k in 1:num.s) {
        for (i in 1:10) {
          grid.lines(x=0:generations, y=sum.hets.indiv[[k]][[i]], gp=gpar(col=k, lwd=2), default.units="native")
        }
      }
      grid.xaxis(gp=gpar(cex=0.8))
      grid.yaxis(gp=gpar(cex=0.8))
    popViewport()
    grid.text("generations", x=0.5, y=-0.18)
    grid.text("avg. sum (over loci) heterozygosity", y=0.5, x=-0.13, rot=90)
  popViewport()
  grid.text(expression(paste("selection: ", Ns, alpha^2)), x=0.825, y=0.92, gp=gpar(cex=0.9))
  pushViewport(viewport(x=0.95, y=0.9, just=c(1,1), height=0.2, width=0.2))
    legend.y <- (1:num.s)/(num.s+1)
    grid.segments(rep(0.1, num.s), legend.y, rep(0.3, num.s), legend.y, gp=gpar(col=1:num.s, lwd=2))
    s <- as.numeric(names(runs.split))
    Nsa2 <- runs[[1]]$par$popsize * s * runs[[1]]$par$effects^2
    grid.text(Nsa2, x=0.4, y=legend.y, just=c(0,0.5), gp=gpar(cex=0.75))
  popViewport()
popViewport()
dev.off()
# END
