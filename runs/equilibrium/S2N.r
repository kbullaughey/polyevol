library(gsl)

N.S2N.est <- function(S, u) 0.5 * S / (u * lambert_W0(1.96895*S/u))
S2N <- function(N, u) 2*N*u*(log(2*N)+0.6775)

u.vec <- c(0.0002, 0.001, 0.01)
N.vec <- 10^(seq(3, 5, length.out=100))

S2N.forward <- lapply(u.vec, function(u) S2N(N.vec, u))
S2N.reverse <- lapply(1:length(u.vec), function(i) sapply(S2N.forward[[i]], function(S) N.S2N.est(S, u.vec[i])))
S2N.rng <- range(unlist(S2N.forward))
N.rng <- range(N.vec)

pdf(file="S2N-check.pdf", height=3, width=4)
par(mar=c(4,4,1,1), mgp=c(2,0.8,0), cex.lab=1, cex.axis=0.75)
plot(log10(N.rng), log10(S2N.rng), type="n", xlab="log10 N", ylab="log10 S")
for (i in 1:length(u.vec)) {
  lines(log10(N.vec), log10(S2N.forward[[i]]))
  lines(log10(S2N.reverse[[i]]), log10(S2N.forward[[i]]), lty=3, col="red", lwd=2)
}
dev.off()

# END

