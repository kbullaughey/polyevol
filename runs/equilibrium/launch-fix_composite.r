N <- c(100,300,1000,3000,10000,30000)
Nu <- c(1,4,16,64,84)
Nsa2 <- c(10,30,100,300)
g <- expand.grid(N=N, Nu=Nu, Nsa2=Nsa2)

for (i in 1:5) {
  cat(paste("qsub run-fix_composite.sh", g$N, g$Nu, g$Nsa2, "\n"), sep="")
}
