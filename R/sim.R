# Constants

VERS <- 2.00
SEED <- 13
NIND <- 1000
TMAX <- 1100  # 110000
BURN <- 100   # 10000
NRUN <- 10  # 100


source("R/functions.R")


b <- 5.0
c <- 0.05
gamma <- 0.1
ar <- double(TMAX+1)
ma <- double(NRUN)
sda <- double(NRUN)
nrep <- integer(NRUN)
set.seed(SEED) # does it work inside functions?
cat("VERS", VERS, "\n")
cat("SEED", SEED, "\n")
cat("NIND", NIND, "\n")
cat("TMAX", TMAX, "\n")
cat("BURN", BURN, "\n")
cat("NRUN", NRUN, "\n")
cat("b", b, "\n")
cat("c", c, "\n")
cat("gamma", gamma, "\n")
cat("#,nrep,ma,sda", "\n")
mnrep <- 0
mma <- 0
msda <- 0
fit <- list(fr = 0, fm = 0) # set initial values???
for(r in 0:(NRUN-1)) {
  ar[1] <- 2.0/(4*c)* runif(1)
  nrep[r+1] <- 0
  for(t in 0:(TMAX-1)) {
    am <- Mutation(ar[t+1],c)
    nm <- 1
    while((nm>0) && (nm<NIND)) {
      fit <- Fitness(b,c,gamma,am,ar[t+1],nm)
      if((fit$fm < 0) || (fit$fr < 0)) {
        cat("negative fitness", r, ", ", t, "\n")
      }
      nm <- Reproduction(nm, fit$fm, fit$fr)
    }
    if(nm == NIND) {
      if(t >= BURN) {
        nrep[r+1] <- nrep[r+1] + 1
      }
      ar[t+2] <- am
    } else {
      ar[t+2] <- ar[t+1]
    }
  }
  ma[r+1] <- 0
  for(i in (BURN+1):TMAX) {
    ma[r+1] <- ma[r+1] + (ar[i+1] / (TMAX-BURN))
  }
  
  sda[r+1] <- 0
  for(i in (BURN+1):TMAX) {
    sda[r+1] <- sda[r+1] + (((ar[i+1]-ma[r+1])^2)/(TMAX-BURN-1))
  }
  sda[r+1] <- sqrt(sda[r+1])
  
  cat(r+2, ", ", nrep[r+1], ", ", ma[r+1], ", ", sda[r+1], "\n")
  mnrep <- mnrep + as.double(nrep[r+1]/NRUN)
  mma <- mma + (ma[r+1] / NRUN)
  msda <- msda + (sda[r+1] / NRUN)
}
cat("mean", " ", mnrep, " ", mma, " ", msda, "\n")
sdnrep <- 0
sdma <- 0
sdsda <- 0
for(r in 0:(NRUN-1)) {
  sdnrep <- sdnrep + (((nrep[r+1]-mnrep)^2)/(NRUN-1.0))
  sdma <- sdma + (((ma[r+1]-mma) ^ 2) / (NRUN-1.0))
  sdsda <- sdsda + (((sda[r+1]-msda) ^ 2) / (NRUN-1.0))
}
sdnrep <- sqrt(sdnrep)
sdma <- sqrt(sdma)
sdsda <- sqrt(sdsda)
cat("sd ", sdnrep, ", ", sdma, ", ", sdsda, "\n")