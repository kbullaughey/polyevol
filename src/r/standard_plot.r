#!/usr/bin/Rscript --verbose

# load arguments
the.args <- commandArgs()
script <- sub("--file=", "", the.args[grep("--file=", the.args)])
run.file <- sub("--data=", "", the.args[grep("--data=", the.args)])
pdf.file <- sub("--pdf=", "", the.args[grep("--pdf=", the.args)])

# check args
stopifnot(length(run.file) == 1)
stopifnot(length(script) == 1)
stopifnot(length(grep("/standard_plot.r$", script)) == 1)

# load library
script.dir <- sub("/standard_plot.r", "", script)
source(paste(script.dir, "/evolveq.R", sep=""))

# assume a name for a pdf file if none is given
if (length(pdf.file) == 0) {
  base.name <- sub("\\.out(\\.gz)*", "", run.file)
  pdf.file <- paste(base.name, ".pdf", sep="")
}

# load the data plot the run
run <- load.simoutput(run.file)
plot.freqs(run, pdf.file=pdf.file, resolution=500)

# END
