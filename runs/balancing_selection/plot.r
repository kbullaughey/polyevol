library(grid)
source("../../src/r/evolveq.R")

pip <- pipe("find out -name 'shift_opt_0_to_1-times_100_1900-s_0.1.*.out'", "r")
files <- scan(pip, what="", sep="\n")
close(pip)

runs <- lapply(files, load.simoutput)

pdf(file="balancing_selection_runs.pdf", height=3.5, width=6)
trash <- lapply(runs, function(run) {
  grid.newpage()
  plot.freqs(run, pdf.file=NULL, color.by="effects")
})
dev.off()


# END
