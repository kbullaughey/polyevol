source("../../src/r/evolveq.R")

run <- "2012_03_13"
load(file="setup-phenotype_effect_grid.rimage")

pip <- pipe("lfs find out/phenotype_effect_grid/ -name 'phenotype_effect_grid*.out'", "r")
files <- scan(pip, what="", sep="\n")
close(pip)

data <- list()
for (i in 1:length(files)) {
  cat(i, "\n")
  # read in the pheno_segsites_compare data
  f <- files[i]
  cat(f, "\n")
  lines <- scan(file=f, what="", sep="\n")
  if (length(lines) < 3) {
    data[[i]] <- NULL
    next
  }
  params <- parse.params(lines)

  # contains a tuple of the sum of phenotypic variances across generations, and the number of generations
  phenovar <- as.numeric(strsplit(grep("phenovarsum:", lines, value=TRUE), " ")[[1]][c(2,4)])
  stopifnot(params$times+1 == phenovar[2])
    
  data[[i]] <- list(params=params, phenovarsum=phenovar[1], count=phenovar[2])
}

save.image(file=paste("phenotype_effect_grid-", run, ".rimage", sep=""))

mx.phenovarsum <- matrix(0, ncol=length(mu), nrow=length(effects))
mx.counts <- matrix(0, ncol=length(mu), nrow=length(effects))

# fill the matrix
for (k in 1:length(data)) {
  x <- data[[k]]
  if (is.null(x)) next
  j <- match(as.character(x$params$mu), as.character(mu))
  i <- match(as.character(x$params$effects), as.character(effects))
  mx.phenovarsum[i,j] <- mx.phenovarsum[i,j] + x$phenovarsum
  mx.counts[i,j] <- mx.counts[i,j] + x$count
}

mx <- mx.phenovarsum / mx.counts

pdf(file="phenotype_effect_contour.pdf", height=4, width=4)
par(mar=c(4,4,3,1), mgp=c(2.2, 0.8, 0), cex.axis=0.75)
contour(2000*mu, effects, t(mx), xlab="mutations/generation", ylab="effect size")
mtext(side=3, line=1.5, "phenotype variance")
dev.off()

# END
