source("../../src/r/evolveq.R")

# read in the set filenames and extract parameters
pip <- pipe("lfs find out/fix_composite/ -D 1 -type d | grep Nsa2", "r")
sets <- scan(pip, what="", sep="\n")
close(pip)
set.param.mx <- matrix(unlist(strsplit(sets, split="-")), ncol=3, byrow=TRUE)
set.param.d <- as.data.frame(apply(sub(".*_", "", set.param.mx), 2, as.numeric))
names(set.param.d) <- sub("_.*$", "", sub(".*/", "", set.param.mx[1,]))
set.param.split <- split(set.param.d, paste(set.param.d$Nu, set.param.d$Nsa2))
sets.split <- split(sets, paste(set.param.d$Nu, set.param.d$Nsa2))

visits.all <- list()
for (i in 1:length(set.param.split)) {
  cat(i, "\n")
  visits.all[[i]] <- lapply(1:length(sets.split[[i]]), function(k) {
    dir.name <- sets.split[[i]][k]
    file.names <- system(paste("ls ", dir.name, "/*.visits.gz", sep=""), intern=TRUE)
    
    # read in the visits data
    visits <- lapply(file.names, function(f) {
      lines <- scan(file=f, what="", sep="\n")
      params <- parse.params(lines)
      lines <- lines[grep("visits: ", lines)]
      if (length(lines) == 1) {
        lines <- sub("^.*visits: ", "", lines)
        v <- as.numeric(strsplit(lines, " ")[[1]])
      } else {
        v <- numeric(0)
      }
      list(visits=v, par=params, file=f)
    })
  })
}

problem.files <- lapply(visits.all, function(x) {
  lapply(x, function(y) {
    lapply(y, function(z) {
      if (length(z$visits) != 2*z$par$popsize-1) {
        return(z$file)
      } else {
        return(NA)
      }
    })
  })
})

save.image(file="fix_composite-visits_distribution_quantiles.rimage")

fbins <- c(0.01, 0.02, 0.03, 0.05, 0.10, 0.5, 0.9)
visits.use <- lapply(visits.all, function(x) {
  sapply(x, function(y) {
    v <- apply(sapply(y, function(z) z$visits), 1, sum)
    freqs <- (1:length(v))/(length(v)+1)
    sapply(fbins, function(b) sum(v[freqs >= b])/sum(v))
  })
})

# sort the plots
params.num <- matrix(as.numeric(unlist(strsplit(names(set.param.split), " "))), ncol=2, byrow=TRUE)
o <- order(params.num[,1], params.num[,2])
visits.use.o <- visits.use[o]
set.param.split.o <- set.param.split[o]

densities.min <- sapply(visits.use.o, function(x) {
  x <- unlist(x)
  min(x[x>0])
})
bins.use <- 1:length(fbins)

pdf(file="fix_composite-frequency_summary.pdf", height=3, width=4)
palette(c("olivedrab", "firebrick", "purple", "tan", "gray30", "orange", "gray60", "dodgerblue", "black"))
par(mar=c(4,4,2.2,2), mgp=c(2,0.8,0), cex.lab=1.0, cex.axis=0.75, xpd=NA, cex.main=0.8)
for (i in 1:length(visits.use.o)) {
  neg.inf <- -11
  N <- set.param.split.o[[i]]$N
  plot(log10(range(N)), c(neg.inf,0), type="n", xlab="N", ylab="log10 density", 
    bty="n", axes=FALSE)
  text(2.5, 1, labels=substitute(N*u == a, list(a=set.param.split.o[[i]]$Nu[1])), cex=0.8) 
  text(3.5, 1, labels=substitute(Ns*alpha^2 == b, list(b=set.param.split.o[[i]]$Nsa2[1])), cex=0.8)
  axis(1, at=log10(N), labels=N)
  vd <- matrix(unlist(visits.use.o[[i]]), ncol=length(N))
  o <- order(N)
  N <- N[o]
  vd <- vd[,o]
  offsets <- apply(vd == 0, 2, function(x) rev(cumsum(rev(x)))*x)*(-neg.inf)*0.03-0.4
  for (k in bins.use) {
    y <- y.orig <- log10(vd[k,])
    off <- offsets[k,]
    y[y == -Inf] <- neg.inf
    stopifnot(sum(y < neg.inf) == 0)
    lines(log10(N), y+off, lty=2, col=k)
    points(log10(N), y+off, pch=20, col=k)
  }
  h <- neg.inf+max(offsets+0.4)
  segments(0.9*log10(N[1]), h, 1.01*log10(max(N)), h)
  at <- seq(neg.inf, 0, by=1)
  at <- at[at > h]
  at <- pretty(at, length=5)
  at <- at[at > h]
  axis(2, at=at, labels=at)

  # calculate positions for the percentile labels
  top <- max(log10(vd[bins.use,ncol(vd)]))
  bottom <- max(neg.inf,min(log10(vd[bins.use,ncol(vd)])))
  total <- (0-neg.inf)
  mid <- (top + bottom)/2
  h <- total*0.26
  n.pos <- seq(mid-h, mid+h, length.out=nrow(vd))
  # plot labels
  n <- paste(fbins[bins.use]*100, "%", sep="")
  text(labels=n, log10(max(N))*1.03, rev(n.pos), cex=0.7, adj=c(0,0.5), col=1:length(bins.use))
}
dev.off()





# END
