
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
qbins <- c(0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.99, 0.999)

sojourn.q.all <- list()
for (i in 1:length(set.param.split)) {
  cat(i, "\n")
  sojourn.q.all[[i]] <- lapply(1:length(sets.split[[i]]), function(k) {
    dir.name <- paste("out/fix_composite/", sets.split[[i]][k], sep="")
    file.names <- system(paste("ls ", dir.name, "/*.sojourns.gz", sep=""), intern=TRUE)
    
    # read in the sojourn data
    soj <- unlist(lapply(file.names, function(f) {
      cat("file:", f, "\n")
      as.numeric(system(paste("zcat", f, "|", "awk '{print $1}'"), intern=TRUE))
    }))
    quantile(soj, probs=qbins)
  })
}
save.image(file="fix_composite-sojourn_distribution_quantiles.rimage")

# log axis
sj.rng <- range(sapply(sojourn.q.all, range))
at <- t(matrix(seq(2, 10, by=2), nrow=4, ncol=10, byrow=TRUE) * 10^(0:3))
at <- c(1,unique(at[1:prod(dim(at))]))
at <- at[at < max(pretty(sj.rng))]
labs <- as.character(at)
labs.sel <- grep("10*", labs)

bins.use <- c(1,3:8)
pdf(file="fix_composite-sojourn_times_distr.pdf", height=3, width=4)
palette(c("olivedrab", "firebrick", "purple", "tan", "gray30", "orange", "gray60", "dodgerblue"))
par(mar=c(4,4,2.2,2), mgp=c(2,0.8,0), cex.lab=1.0, cex.axis=0.75, xpd=NA, cex.main=0.8)
for (i in 1:length(sojourn.q.all)) {
  N <- set.param.split[[i]]$N
  plot(log10(range(N)), log10(sj.rng), type="n", xlab="N", ylab="sojourn time", 
    bty="n", axes=FALSE)
  text(2.5, 3.5, labels=substitute(N*u == a, list(a=set.param.split[[i]]$Nu[1])), cex=0.8) 
  text(3.5, 3.5, labels=substitute(Ns*alpha^2 == b, list(b=set.param.split[[i]]$Nsa2[1])), cex=0.8)
  axis(1, at=log10(N), labels=N)
  axis(2, at=log10(at), labels=FALSE)
  text(1.77, log10(at[labs.sel]), labels=labs[labs.sel], cex=0.6, adj=c(1,0.5))
  sjq <- matrix(unlist(sojourn.q.all[[i]]), ncol=length(N))
  o <- order(N)
  N <- N[o]
  sjq <- sjq[,o]
  for (k in bins.use) {
    lines(log10(N), log10(sjq[k,]), lty=2, col=k)
    points(log10(N), log10(sjq[k,]), pch=20, col=k)
  }
  n <- names(sojourn.q.all[[1]][[1]])
  text(labels=n[bins.use], log10(max(N))*1.03, log10(sjq[bins.use,ncol(sjq)]), cex=0.7, adj=c(0,0.5))
}
dev.off()





# END
