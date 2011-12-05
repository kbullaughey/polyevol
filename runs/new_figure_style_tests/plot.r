sets <- c("balancing-infinite", "balancing-finite", "no_balancing-infinite", "no_balancing-finite")
sets <- c("no_balancing-infinite", "no_balancing-finite")

set.runs <- lapply(sets, function(s) {
  # load the runs
  pip <- pipe(paste("ls ", s, ".*.out", sep=""), "r")
  files <- scan(file=pip, what="", sep="\n")
  close(pip)
  runs <- lapply(files, load.simoutput)

  # plot the runs
  pdf(file=paste(s, ".pdf", sep=""), height=5.5, width=6)
  for (i in 1:length(runs)) {
    if (i > 1) grid.newpage()
    plot.freqs(runs[[i]], pdf.file=NULL, resolution=1000)
  }
  dev.off()
})

# END
