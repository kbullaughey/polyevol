# Load empirical data
pip <- pipe("lfs find out -name 'summary' | grep 'phenotype-N'", "r")
summary.files <- scan(pip, what="", sep="\n")
close(pip)

N <- as.numeric(sub("^.*N_([0-9]*)-.*$", "\\1", summary.files))
s <- 1/as.numeric(sub("^.*s_([.0-9]*)/.*$", "\\1", summary.files))
combinations <- data.frame(N=N, s=s)

summaries <- lapply(summary.files, function(f) {
  tbl <- read.table(f, header=FALSE)
  names(tbl) <- c("mu", "pheno", "count")
  tbl <- as.data.frame(t(sapply(split(tbl, tbl$mu), function(x) c(mu=x$mu[1], pheno=sum(x$pheno)/sum(x$count)))))
  rownames(tbl) <- NULL
  return(tbl)
})
empirical.end.points <- sapply(summaries, function(x) rev(x$pheno)[1])

f.diploid <- function(v, U, N, s) 8*N*U * ((v+s)/(2*N) - (sqrt(1+v/s*v/(2*v+s))*(exp(2*N/(v+s)/sqrt(1+v/s*v/(2*v+s))) - 1))^(-1)) - v

u <- seq(0, 0.5, length.out=100)

v <- lapply(1:nrow(combinations), function(i) {
  sapply(u, function(ui) uniroot(f.diploid, c(0,10^6), ui, combinations$N[i], combinations$s[i])$root)
})
true.end.points <- sapply(v, function(x) rev(x)[1])

# reorder by end points
o <- order(empirical.end.points)
v <- v[o]
true.end.points <- true.end.points[o]
summaries <- summaries[o]
combinations <- combinations[o,]
empirical.end.points <- empirical.end.points[o]

v.rng <- range(pretty(unlist(v)))

pdf(file="phenotypic_variance-mean_field_compare.pdf", height=5, width=6.5)
par(mar=c(4,4,1,5), mgp=c(2,0.8,0), cex.axis=0.75, xpd=NA)
palette(c("black", "orange", "firebrick", "olivedrab", "gray40", "dodgerblue", "tan", "darkseagreen", "darksalmon", "purple"))
plot(range(u), c(0,800), type="n", xlab="mutation rate (U)", ylab="phenotypic variance (v)", bty="n")
trash <- lapply(1:nrow(combinations), function(i) {
  lines(u, v[[i]], lwd=1, col=i)
  points(summaries[[i]]$mu, summaries[[i]]$pheno, pch=20, col=i)
})
label <- paste("N: ", combinations$N, ", s: ", combinations$s, sep="")

# Here I calculate positions for the labels using an optimization algorithm 
# to determine the beset label placement
objective <- function(x, p, h) sum(abs(cumsum(x+h) - h/2 - p))
initial <- rep(0, nrow(combinations))
h <- v.rng[2]*0.03
res <- optim(initial, objective, gr=NULL, empirical.end.points, h, 
  method="L-BFGS-B", lower=0)
n.pos <- cumsum(res$par + h) - h/2
segments(0.41, empirical.end.points, 0.43, n.pos, lty=3, col=the.colors)
text(0.44, n.pos, labels=label, cex=0.6, col=the.colors, adj=c(0,0.5))
dev.off()
