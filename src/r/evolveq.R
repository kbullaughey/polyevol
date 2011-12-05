library(grid)

# reduce a long vector to an evenly spaced subset
reduce.index <- function(x, len=100) unique(floor(seq(1, length(x), length.out=len)))
reduce <- function(x, len=100) x[reduce.index(x,len=len)]

smush <- function(x, sep=",", digits=2) paste(signif(x, digits=digits), collapse=sep)

# mix up a vector
shuffle <- function(x) x[sample(1:length(x), length(x), replace=FALSE)]

# don't transform
pass.thru <- function(x, N) x

# log10 in both directions away from 0.5 mapped onto to interval (1-1/2N, 1/2N) -> (0,1)
bilog <- function(x, N) {
  res <- numeric(length(x))
  boundary <- x < 1/(2*N) | x > (1-1/(2*N))
  high <- x > 0.5 & !boundary
  low <- x <= 0.5 & !boundary
  res[low] <- (log10(x[low]) - log10(1/(2*N))) / (2*(log10(0.5) - log10(1/(2*N))))
  res[high] <- 1 - (log10(1-x[high]) - log10(1/(2*N))) / (2*(log10(0.5) - log10(1/(2*N))))
  res[boundary] <- NA
  return(res)
}

# Extract and use the frequency data from the simulation output to built trajectories
# for each locus
parse.trajectories.from.frequency.data <- function(x) {
  freq.lines <- x[grep("freqs:", x)]
  # extract the generation from each line
  gen <- as.numeric(sub("^.*gen: ([0-9]*) .*$", "\\1", freq.lines))
  # split the lines on spaces, now a list of character vectors of the form "id:freq"
  x <- strsplit(sub("^.*freqs: ", "", freq.lines), split=" ")
  # extract out the IDs and frequencies
  locus.ids <- lapply(x, function(z) as.numeric(sapply(strsplit(z, split=":"), function(y) y[1])))
  locus.freqs <- lapply(x, function(z) as.numeric(sapply(strsplit(z, split=":"), function(y) y[2])))
  # the number of loci each generation, which can vary in the infinite sites model
  num.loci <- sapply(locus.ids, length)
  # make sure our stat structures are the same size
  stopifnot(sum(sapply(locus.freqs, length) != num.loci) == 0)
  # put all this into a data frame with three columns: gen, locus, freq
  df <- data.frame(gen=rep(gen, num.loci), locus=unlist(locus.ids), freq=unlist(locus.freqs))
  # also split things by loucs:
  df.list <- split(df[,-match("locus", names(df))], df$locus)
  # return both
  return(list(df=df, lst=df.list))
}

# Turn a deliminated bunch of lines into a matrix
parse.to.matrix <- function(x, sep=" ", transform=as.numeric) {
  spl <- lapply(strsplit(x, split=sep), transform)
  # require all lines to have the same number of entries
  stopifnot(length(unique(sapply(spl, length))) == 1)
  matrix(unlist(spl), ncol=length(spl[[1]]), byrow=TRUE)
}

# load an output file from the simulation into a friendly package
load.simoutput <- function(file.name) {
  file.lines <- scan(file=file.name, what="", sep="\n")

  # get the raw command
  cmd <- sub("^cmd: ", "", file.lines[grep("^cmd: ", file.lines)])

  # use the params line to get all configuration
  params <- sub("^params: ", "", file.lines[grep("^params: ", file.lines)])
  stopifnot(length(params) == 1)
  params <- sub(",$", "", strsplit(params, split=" ")[[1]])
  # make a list out of the parameters, which are formatted as R assignments
  params <- eval(parse(text=paste("list(", paste(params, collapse=","), ")")))

  # extract the data from the file in subsets based on line prefixes
  traj <- parse.trajectories.from.frequency.data(file.lines)

  gen.and.pheno <- sub("^.*gen: ([0-9]+) pheno: ", "\\1 ", file.lines[grep("pheno:", file.lines)])
  phenodata <- as.data.frame(parse.to.matrix(gen.and.pheno))
  stopifnot(ncol(phenodata) == 3)
  names(phenodata) <- c("gen", "mean", "var")

  # extract the effect sizes from the 'site' output lines
  site.lines <- file.lines[grep("^site:", file.lines)]
  site.id <- sub("^.* id: ([0-9]+) .*$", "\\1", site.lines)
  site.effect <- as.numeric(sub("^.*effect: (-?[0-9]).*$", "\\1", site.lines, perl=TRUE))
  names(site.effect) <- site.id

  list(freq=traj$lst, pheno=phenodata, par=params, effects=site.effect, cmd=cmd)
}

# pass pdf.file=NULL in order to use an open device
# color.by can take two values: 
#   "all"       - each locus colored a different color 
#   "effects"   - each effect size colored a different color
# pass.thru is used to transform the frequencies onto a different scale for plotting it can 
#   take either a string naming a function or the function itself. Supported functions are:
#   - pass.thru : don't do anything, just use the frequencies directly
#   - bilog : log10 in both directions away from 0.5
#   the function given or named must take two arguments, a vector to transform and N
plot.freqs <- function(run, pdf.file="out.pdf", resolution=500, color.by="all", freq.tr="pass.thru") {
  nloci <- length(run$effects)

  # plotting
  if (!is.null(pdf.file))
    pdf(file=pdf.file, height=5.5, width=6)

  # determine coloring
  if (color.by == "all") {
    color.id <- 1:nloci %% 25
  } else if (color.by == "effects") {
    color.id <- as.numeric(factor(abs(run$effects)))
  } else {
    stop("invalid color.by")
  }
  names(color.id) <- names(run$effects)
  num.colors <- length(unique(color.id))
  if (num.colors <= 5) {
    palette(c("purple", "firebrick", "tan", "steelblue", "orange"))
  } else {
    palette(shuffle(rainbow(num.colors)))
  }

  x <- strsplit(run$cmd, split=" (?=-)", perl=TRUE)[[1]]
  title.lines <- sapply(split(x, floor(cumsum(nchar(x)) / 70)), paste, collapse=" ")
  title.lines <- paste(
    c("", rep("   ", length(title.lines)-1)), 
    title.lines, 
    c(rep(" \\", length(title.lines)-1), ""), sep="")
  title.y <- 0.80-((1:length(title.lines))-1)/4

  # plot the frequencies
  pushViewport(viewport(width=0.75))
    pushViewport(viewport(layout=grid.layout(nrow=3, ncol=1, heights=c(0.1, 0.45, 0.45))))
      pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))
        grid.text(title.lines, x=0.03, y=title.y, gp=gpar(cex=0.7, fontfamily="mono"), just=c(0,1))
      popViewport()

      pushViewport(viewport(layout.pos.row=2, layout.pos.col=1))
        plot.freqs.aux(run, color.id, resolution=resolution, freq.tr="pass.thru")
      popViewport()

      pushViewport(viewport(layout.pos.row=3, layout.pos.col=1))
        y.at <- signif(1/10^seq(0,floor(log10(run$par$popsize)), by=0.5)/2, digits=2)
        y.at <- unique(c(rev(1-y.at), y.at))
        plot.freqs.aux(run, color.id, resolution=resolution, freq.tr="bilog", at=y.at)
      popViewport()
    popViewport()
  popViewport()

  # close up shop
  if (!is.null(pdf.file))
    dev.off()
}

plot.freqs.aux <- function(run, color.id, resolution=500, freq.tr="pass.thru", at=seq(0,1,by=0.2)) {
  ngen <- tail(run$par$times,1)
  pushViewport(viewport(y=1, just=c(0.5,1), height=0.8, xscale=c(0,ngen), yscale=0:1))
    grid.segments(run$par$times, 0, run$par$times, 1, default.units="native", gp=gpar(col="gray", lty=3, lwd=2))
    trash <- lapply(1:length(run$freq), function(i) {
      id <- names(run$freq)[i]
      g <- run$freq[[i]]$gen
      if (run$effects[i] < 0) {
        f <- 1-run$freq[[i]]$freq
      } else {
        f <- run$freq[[i]]$freq
      }
      f <- do.call(freq.tr, list(f, run$par$popsize))
      grid.lines(g, f, gp=gpar(col=color.id[[id]]), default.units="native") 
    })
    grid.xaxis(gp=gpar(cex=0.7))
    y.at <- do.call(freq.tr, list(at, run$par$popsize))
    grid.yaxis(label=at, at=y.at, gp=gpar(cex=0.7))
    grid.text("Generations", x=0.5, y=-0.18, gp=gpar(cex=0.9))
    grid.text("Frequencies", y=0.5, x=-0.12, rot=90, gp=gpar(cex=0.9))
  popViewport()

  # plot the phenotype
  y <- run$pheno$mean
  g <- run$pheno$gen
  pushViewport(viewport(y=1, just=c(0.5,1), height=0.8, xscale=c(0,ngen), yscale=range(c(y, run$par$opts))))
    idx <- reduce.index(y, len=resolution)
    grid.lines(g[idx], y[idx], gp=gpar(col=rgb(0,0,0,0.3), lwd=3), default.units="native") 
    grid.yaxis(gp=gpar(cex=0.75), main=FALSE)
    grid.text("Phenotype", y=0.5, x=1.1, rot=90, gp=gpar(cex=0.9))
    grid.rect()
  popViewport()
}

# END
