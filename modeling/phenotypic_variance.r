# Load empirical data
pip <- pipe("lfs find ../runs/equilibrium/out -name 'summary' | grep 'phenotype-N'", "r")
summary.files <- scan(pip, what="", sep="\n")
close(pip)

f.haploid <- function(v, U, N, s) 2*N*U * ((v+s)/N - (sqrt(1+v/s*v/(2*v+s))*(exp(N/(v+s)*sqrt(1 - v^2/(v+s)^2)) - 1))^(-1)) - v
f.diploid <- function(v, U, N, s) 8*N*U * ((v+s)/(2*N) - (sqrt(1+v/s*v/(2*v+s))*(exp(2*N/(v+s)/sqrt(1+v/s*v/(2*v+s))) - 1))^(-1)) - v

N <- c(100, 1000, 3000)
s <- c(10, 30, 100, 200)
combinations <- expand.grid(N=N, s=s)

u <- seq(0, 1, length.out=100)

v <- lapply(1:nrow(combinations), function(i) {
  sapply(u, function(ui) uniroot(f.diploid, c(0,10^6), ui, combinations$N[i], combinations$s[i])$root)
})

v.rng <- range(pretty(unlist(v)))

pdf(file="phenotypic_variance-mean_field_approx-diploid.pdf", height=5, width=6.5)
par(mar=c(4,4,1,4), mgp=c(2,0.8,0), cex.axis=0.75, xpd=NA)
palette(c("black", "orange", "gray", "firebrick", "olivedrab", "dodgerblue", "tan", "darkseagreen", "darksalmon", "purple"))
plot(range(u), v.rng, type="n", xlab="mutation rate (U)", ylab="phenotypic variance (v)", bty="n")
trash <- lapply(1:nrow(combinations), function(i) {
  lines(u, v[[i]], lwd=1, col=i)
})
label <- paste("N: ", combinations$N, ", s: ", combinations$s, sep="")
true.end.points <- sapply(v, function(x) rev(x)[1])

# reorder by end points
o <- order(true.end.points)
v <- v[o]
true.end.points <- true.end.points[o]
the.colors <- (1:nrow(combinations))[o]

# Here I calculate positions for the labels using an optimization algorithm 
# to determine the beset label placement
objective <- function(x, p, h) sum(abs(cumsum(x+h) - h/2 - p))
initial <- rep(0, nrow(combinations))
h <- v.rng[2]*0.03
res <- optim(initial, objective, gr=NULL, true.end.points, h, 
  method="L-BFGS-B", lower=0)
n.pos <- cumsum(res$par + h) - h/2
segments(1.01, true.end.points, 1.05, n.pos, lty=3, col=the.colors)
text(1.05, n.pos, labels=label, cex=0.6, col=the.colors, adj=c(0,0.5))
dev.off()
