source("../../src/r/evolveq.R")

the.args <- commandArgs()
file.arg <- sub("--file=", "", the.args[grep("--file=", the.args)])
file.freq <- sub(".out$", ".freq_distr", file.arg)

run <- load.simoutput(file.arg)

distr <- numeric(run$par$popsize*2)
tbl <- table(unlist(lapply(run$freq, function(x) x$freq)))
distr[as.numeric(names(tbl)) * run$par$popsize*2] <- as.numeric(tbl)
write(distr, file=file.freq, ncolumns=10)

# END
