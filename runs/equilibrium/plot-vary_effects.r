source("../../src/r/evolveq.R")

run <- "2012_03_06"

# read in the set filenames and extract parameters
pip <- pipe("lfs find out/vary_effects/ -D 1 -type d | grep mu", "r")
sets <- scan(pip, what="", sep="\n")
close(pip)

summaries <- list()
for (i in 1:length(sets)) {
  cat(i, "\n")
  dir.name <- sets[i]
  file.names <- system(paste("lfs find ", dir.name, " -name vary_effects*.out.gz", sep=""), intern=TRUE)
    
  # read in the fixations data
  summaries[[i]] <- lapply(file.names, function(f) {
    cat(f, "\n")
    lines <- scan(file=f, what="", sep="\n")
    params <- parse.params(lines)
    pheno.var <- as.numeric(sapply(strsplit(grep("pheno:", lines, value=TRUE), " "), function(x) x[5]))
#    visits <- as.numeric(strsplit(sub("visits: ", "", grep("visits:", lines, value=TRUE)), split=" ")[[1]])
#    stopifnot(length(visits) == params$popsize*2-1)
    list(params=params, pheno.var=pheno.var)
  })
}

save.image(file=paste("vary_effects-", run, ".rimage", sep=""))


# END
