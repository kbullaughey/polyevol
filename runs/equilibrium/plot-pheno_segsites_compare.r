source("../../src/r/evolveq.R")

run <- "2012_03_07"

# read in the set filenames
pip <- pipe("lfs find out/pheno_segsites_compare/ -D 1 -type d | grep 'u-'", "r")
sets <- scan(pip, what="", sep="\n")
close(pip)

data <- list()
for (i in 1:length(sets)) {
  cat(i, "\n")
  dir.name <- sets[i]
  file.names <- system(paste("lfs find ", dir.name, " -name pheno_segsites_compare*.out", sep=""), intern=TRUE)
    
  # read in the pheno_segsites_compare data
  data[[i]] <- lapply(file.names, function(f) {
    cat(f, "\n")
    lines <- scan(file=f, what="", sep="\n")
    params <- parse.params(lines)
    # contains a tuple of the sum of neutral segregating sites, and the number of generations
    neutral <- as.numeric(strsplit(grep("neutral:", lines, value=TRUE), " ")[[1]][-1])
    stopifnot(params$times+1 == neutral[2])
    
    # contains a tuple of the sum of phenotypic variances across generations, and the number of generations
    phenovar <- as.numeric(strsplit(grep("phenovar:", lines, value=TRUE), " ")[[1]][-1])
    stopifnot(params$times+1 == phenovar[2])
    
    list(params=params, neutral=neutral, phenovar=phenovar)
  })
}

save.image(file=paste("pheno_segsites_compare-", run, ".rimage", sep=""))

# collate the data
data.col <- as.data.frame(t(sapply(data, function(z) {
  neutral.combined <- apply(sapply(z, function(x) x$neutral), 1, sum)
  phenovar.combined <- apply(sapply(z, function(x) x$phenovar), 1, sum)
  mu <- unique(sapply(z, function(x) x$params$mu))
  stopifnot(length(mu) == 1)
  c(mu=mu, neutral=neutral.combined[1]/neutral.combined[2], phenovar=phenovar.combined[1]/phenovar.combined[2])
})))
data.col <- data.col[order(data.col$mu),]

# estimate N
N.hat <- sapply(1:nrow(data.col), function(i) N.S2N.est(data.col$neutral[i], data.col$mu[i]))

# use non-neutral mutations per generation (half the mutations on average have effect 1)
mu <- data.col$mu * 2000 / 2

mu.rng <- range(pretty(mu))
phenovar.rng <- range(pretty(data.col$phenovar))
N.rng <- range(pretty(c(0, N.hat)))
big.mar <- 0.22
small.mar <- 0.04

pdf(file="pheno_segsites_compare.pdf", height=3, width=5)
# plot phenotype variance as function of mutation rate
pushViewport(viewport(x=1-big.mar, y=1-small.mar, height=1-big.mar-small.mar, width=1-2*big.mar, 
  just=c(1,1), xscale=mu.rng, yscale=phenovar.rng))
grid.lines(x=mu, y=data.col$phenovar, default.units="native")
grid.xaxis(gp=gpar(cex=0.8))
grid.yaxis(gp=gpar(cex=0.8))
grid.text("phenotypic variance", y=0.5, x=-0.22, rot=90)
grid.text("non-neutral mutations/gen", x=0.5, y=-0.22)
popViewport()
# plot phenotype variance as function of mutation rate
pushViewport(viewport(x=1-big.mar, y=1-small.mar, height=1-big.mar-small.mar, width=1-2*big.mar, 
  just=c(1,1), xscale=mu.rng, yscale=N.rng))
grid.lines(x=mu, y=N.hat, default.units="native", gp=gpar(col="brown"))
grid.yaxis(gp=gpar(cex=0.8, col="brown"), main=FALSE)
grid.text(expression(hat(N)), y=0.5, x=1.22, gp=gpar(col="brown"))
popViewport()
dev.off()

# END
