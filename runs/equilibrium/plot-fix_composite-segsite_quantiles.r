
# read in the set filenames and extract parameters
pip <- pipe("lfs find out/fix_composite/ -D 1 -type d | grep Nsa2", "r")
sets <- scan(pip, what="", sep="\n")
close(pip)
set.param.mx <- matrix(unlist(strsplit(sets, split="-")), ncol=3, byrow=TRUE)
set.param.d <- as.data.frame(apply(sub(".*_", "", set.param.mx), 2, as.numeric))
names(set.param.d) <- sub("_.*$", "", set.param.mx[1,])
Nu.choices <- unique(set.param.d$Nu)
Nsa2.choices <- unique(set.param.d$Nsa2)
set.param.split <- split(set.param.d, paste(set.param.d$Nu, set.param.d$Nsa2))
sets.split <- split(sets, paste(set.param.d$Nu, set.param.d$Nsa2))

#params

segsites.q.all <- list()
for (i in 1:length(set.param.split)) {
  cat(i, "\n")
  segsites.q.all[[i]] <- lapply(1:length(sets.split[[i]]), function(k) {
    dir.name <- sets.split[[i]][k]
    file.names <- system(paste("ls ", dir.name, "/*.segsites.gz", sep=""), intern=TRUE)
    
    # read in the segsite data
    segs <- lapply(file.names, function(f) {
      cat("file:", f, "\n")
      scan(file=f, what=0, sep="\n")
    })
    list(segs=segs, counts=sapply(segs, length), files=file.names)
  })
}
save.image(file="fix_composite-segsite_distribution_quantiles.rimage")

problem.files <- lapply(segsites.q.all, function(x) {
  lapply(x, function(y) {
    y$files[y$counts < 100001]
  })
})

qbins <- c(0.001, 0.05, 0.25, 0.5, 0.75, 0.95, 0.999)
segsites.q.use <- lapply(segsites.q.all, function(x) {
  lapply(x, function(y) {
    quantile(unlist(y$segs[y$counts == 100001]), prob=qbins)
  })
})

# sort the plots
params.num <- matrix(as.numeric(unlist(strsplit(names(set.param.split), " "))), ncol=2, byrow=TRUE)
o <- order(params.num[,1], params.num[,2])
segsites.q.use.o <- segsites.q.use[o]
set.param.split.o <- set.param.split[o]

sj.rngs <- t(sapply(segsites.q.use.o, function(x) range(unlist(x))))
bins.use <- 1:length(qbins)
pdf(file="fix_composite-segsites_distr.pdf", height=3, width=4)
palette(c("olivedrab", "firebrick", "purple", "tan", "gray30", "orange", "gray60", "dodgerblue", "black"))
par(mar=c(4,4,2.2,2), mgp=c(2,0.8,0), cex.lab=1.0, cex.axis=0.75, xpd=NA, cex.main=0.8)
for (i in 1:length(segsites.q.use.o)) {
  N <- set.param.split.o[[i]][[1]]
  plot(log10(range(N)), c(0,sj.rngs[i,2]), type="n", xlab="N", ylab="segregating sites", 
    bty="n", axes=FALSE)
  h <- sj.rngs[i,2] + (sj.rngs[i,2] - sj.rngs[i,1])*0.1
  text(2.5, h, labels=substitute(N*u == a, list(a=set.param.split.o[[i]]$Nu[1])), cex=0.8) 
  text(3.5, h, labels=substitute(Ns*alpha^2 == b, list(b=set.param.split.o[[i]]$Nsa2[1])), cex=0.8)
  axis(1, at=log10(N), labels=N)
  axis(2)
  sjq <- matrix(unlist(segsites.q.use.o[[i]]), ncol=length(N))
  o <- order(N)
  N <- N[o]
  sjq <- sjq[,o]
  for (k in bins.use) {
    lines(log10(N), sjq[k,], lty=2, col=k)
    points(log10(N), sjq[k,], pch=20, col=k)
  }

  # calculate positions for the percentile labels
  top <- max(sjq[bins.use,ncol(sjq)])
  bottom <- min(sjq[bins.use,ncol(sjq)])
  if ((top-bottom)/(sj.rngs[i,2] - sj.rngs[i,1]) < 0.45) {
    mid <- (top + bottom)/2
    h <- (sj.rngs[i,2] - sj.rngs[i,1])*0.26
    n.pos <- seq(mid-h, mid+h, length.out=nrow(sjq))
  } else {
    n.pos <- sjq[bins.use,ncol(sjq)]
  }
  # plot labels
  n <- names(segsites.q.use.o[[1]][[1]])
  text(labels=n[bins.use], log10(max(N))*1.03, n.pos, cex=0.7, adj=c(0,0.5), col=1:length(bins.use))
}
dev.off()





# END
