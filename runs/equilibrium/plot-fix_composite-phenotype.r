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

phenotype.all <- list()
for (i in 1:length(set.param.split)) {
  cat(i, "\n")
  phenotype.all[[i]] <- lapply(1:length(sets.split[[i]]), function(k) {
    dir.name <- sets.split[[i]][k]
    file.names <- system(paste("ls ", dir.name, "/*.pheno.gz", sep=""), intern=TRUE)
    
    # read in the penotype data
    pheno <- lapply(file.names, function(f) {
      cat(f, "\n")
      tbl <- read.table(file=f, header=FALSE)
      names(tbl) <- c("mean", "var")
      return(tbl)
    })
  })
}

problem.files <- lapply(phenotype.all, function(x) {
  sapply(x, function(y) {
    sapply(y, nrow)
  })
})

save.image(file="fix_composite-pheno.rimage")

qbins <- c(0.001, 0.05, 0.25, 0.5, 0.75, 0.95, 0.999)
phenotype.use <- lapply(phenotype.all, function(x) {
  do.call(rbind, lapply(x, function(y) {
    p <- do.call(rbind, y)$var
    quantile(p, prob=qbins)
  }))
})

# sort the plots
params.num <- matrix(as.numeric(unlist(strsplit(names(set.param.split), " "))), ncol=2, byrow=TRUE)
o <- order(params.num[,1], params.num[,2])
phenotype.use.o <- phenotype.use[o]
set.param.split.o <- set.param.split[o]

sj.rngs <- t(sapply(phenotype.use.o, function(x) range(unlist(x))))
bins.use <- 1:length(qbins)
pdf(file="fix_composite-phenotype_variance_quantiles.pdf", height=3, width=4)
palette(c("olivedrab", "firebrick", "purple", "tan", "gray30", "orange", "gray60", "dodgerblue", "black"))
par(mar=c(4,4,2.2,2), mgp=c(2,0.8,0), cex.lab=1.0, cex.axis=0.75, xpd=NA, cex.main=0.8)
for (i in 1:length(phenotype.use.o)) {
  N <- set.param.split.o[[i]]$N
  plot(log10(range(N)), c(0,sj.rngs[i,2]), type="n", xlab="N", ylab="phenotype variance", 
    bty="n", axes=FALSE)
  h <- sj.rngs[i,2] + (sj.rngs[i,2] - sj.rngs[i,1])*0.1
  text(2.5, h, labels=substitute(N*u == a, list(a=set.param.split.o[[i]]$Nu[1])), cex=0.8) 
  text(3.5, h, labels=substitute(Ns*alpha^2 == b, list(b=set.param.split.o[[i]]$Nsa2[1])), cex=0.8)
  axis(1, at=log10(N), labels=N)
  axis(2)
  sjq <- t(phenotype.use.o[[i]])
  o <- order(N)
  N <- N[o]
  sjq <- sjq[,o]
  for (k in bins.use) {
    lines(log10(N), sjq[k,], lty=2, col=k)
    points(log10(N), sjq[k,], pch=20, col=k)
  }

  # calculate positions for the percentile labels
  # I use an optimization algorithm to determine the beset label placement
  objective <- function(x, p, h) sum(abs(cumsum(x+h) - h/2 - p))
  initial <- rep(0, nrow(sjq))
  end.points <- sjq[bins.use,ncol(sjq)]
  h <- sj.rngs[i,2]*0.5/nrow(sjq)
  res <- optim(initial, objective, gr=NULL, 
    end.points, h, method="L-BFGS-B", lower=0)
  n.pos <- cumsum(res$par + h) - h/2

  # plot labels
  n <- colnames(phenotype.use.o[[1]])
  text(labels=n[bins.use], log10(max(N))*1.03, n.pos, cex=0.7, adj=c(0,0.5), col=1:length(bins.use))
}
dev.off()





# END



# END
