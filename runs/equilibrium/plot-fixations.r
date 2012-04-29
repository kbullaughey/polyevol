source("../../src/r/evolveq.R")

run <- "2012_03_01"

# read in the set filenames and extract parameters
pip <- pipe("lfs find out/fixations/ -D 1 -type d | grep effects", "r")
sets <- scan(pip, what="", sep="\n")
close(pip)

# extract parameters
set.mu <- as.numeric(gsub("^out/fixations/u-(.*)--effects_1$", "\\1", sets))

# sort by mutation rate
o <- order(set.mu)
set.mu <- set.mu[o]
sets <- sets[o]

fixations.all <- list()
for (i in 12:length(sets)) {
  cat(i, "\n")
  dir.name <- sets[i]
  file.names <- system(paste("lfs find ", dir.name, " -name fixations*.out.gz", sep=""), intern=TRUE)
    
  # read in the fixations data
  fixations.all[[i]] <- lapply(file.names, function(f) {
    cat(f, "\n")
    lines <- scan(file=f, what="", sep="\n")
    params <- parse.params(lines)
    tuples <- strsplit(sub("gen: [0-9]* fixations: ?", "", grep("fixations:", lines, value=TRUE)), split=" ")
    diff <- sapply(tuples, function(x) {
      sum(apply(matrix(as.numeric(unlist(strsplit(x, ","))), ncol=2, byrow=TRUE), 1, prod))
    })
    list(params=params, diff=diff)
  })
}

save.image(file=paste("fixations-", run, ".rimage", sep=""))

# compute tables tallying the frequency of various fixation differences
fixations.combined <-lapply(fixations.all, function(y) do.call(c, lapply(y, function(x) x$diff)))
fixations.tables <- lapply(fixations.combined, table)
x.range <- range(unlist(sapply(fixations.tables, function(x) as.numeric(names(x)))))
all <- unlist(sapply(fixations.tables, function(x) as.numeric(x)))
y.max <- max(all)
y.nonzero.min <- min(all[all > 0])
num.sets <- length(sets)

# log axis
high <- ceiling(log10(y.max))
at <- t(matrix(seq(2, 10, by=2), nrow=high, ncol=5, byrow=TRUE) * 10^(1:high - 1))
at <- c(1,unique(at[1:prod(dim(at))]))
labs <- as.character(at)
labs.sel <- grep("10*", labs)

# plot
pdf(file=paste("fixations-", run, ".pdf", sep=""), height=3, width=4.5)
palette(c("black", "olivedrab", "firebrick", "purple", "tan", "gray30", "orange", "gray60", "dodgerblue", "pink", "lightblue"))
par(mar=c(4,4,2.2,4), mgp=c(2,0.7,0), cex.lab=1.0, cex.axis=0.65, xpd=NA, cex.main=0.8)
plot(range(pretty(x.range)), log10(c(y.nonzero.min, y.max))+c(-1,0), type="n", ylab="", xlab="fixation state", 
  bty="n", axes=FALSE)
axis(1)
axis(2, at=log10(at), labels=FALSE)
mtext(side=2, line=2.6, "log10 frequency")
text(x=x.range[1]-0.11*(x.range[2]-x.range[1]), y=log10(at[labs.sel]), labels=labs[labs.sel], adj=c(1, 0.5), cex=0.6)
for (i in 1:num.sets) {
  x <- as.numeric(names(fixations.tables[[i]]))
  y <- as.numeric(fixations.tables[[i]])
  lines(x, log10(y), col=i)
  points(x, log10(y), col=i, pch=21, bg=i, cex=0.6)
}
# legend
legend((x.range[2]-x.range[1])*0.17+x.range[2], high, bty="n", title="mutations/gen", yjust=1, xjust=0, x.intersp=2.2,
  legend=2*1000*set.mu, col=1:num.sets, lwd=1, cex=0.6, adj=c(1, 0.5))
dev.off()





# END
