source("../r/evolveq.R")

run <- load.simoutput("run.out")
plot.freqs(run, pdf.file="run.pdf", color.by="effects", resolution=200)

# END
