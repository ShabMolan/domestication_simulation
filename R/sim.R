VERS <- 2.00
SEED <- 13
NIND <- 1000
TMAX <- 1100  # 110000
BURN <- 100   # 10000
NRUN <- 10  # 100

Mutation <- function(ar, c) { # why c and what is ar
  return(2*1.0/(4*c) * runif(1))
}


Fitness <- function(b, c, gamma, am, a, k) { # removed the *
  if(b<a*(a+3*am)/((a+am)*(a-am))){
    fm <- 1+gamma*(HGD(0,NIND-1,k-1,2)*am*(a+3*am)/(3*(a+am)*(a+am))
                   +HGD(1,NIND-1,k-1,2)*am*(5*a+3*am)/(6*(a+am)*(a+am))
                   +HGD(2,NIND-1,k-1,2)/3.0-c*am)
    
    fr <- 1+gamma*(HGD(0,NIND-1,k,2)/3.0
                   +HGD(1,NIND-1,k,2)*a*(3*a+5*am)/(6*(a+am)*(a+am))
                   +HGD(2,NIND-1,k,2)*a*(3*a+am)/(3*(a+am)*(a+am))-c*a)
    
  }else if((b>a*(a+3*am)/((a+am)*(a-am))) && (b<(3*a+am)/(2*(a-am)))){
    fm <- 1+gamma*(HGD(0,NIND-1,k-1,2)*(2*am/(3*(a+am))+b*am/(3*(b*(a+am)+a)))
                   +HGD(1,NIND-1,k-1,2)*am*(5*a+3*am)/(6*(a+am)*(a+am))
                   +HGD(2,NIND-1,k-1,2)/3.0-c*am)
    
    fr <- 1+gamma*(HGD(0,NIND-1,k,2)/3.0
                   +HGD(1,NIND-1,k,2)*(a/(3*(a+am))+(1+b)*a/(6*(b*(a+am)+a)))
                   +HGD(2,NIND-1,k,2)*a*(3*a+am)/(3*(a+am)*(a+am))-c*a)
    
  }else if(2*b*(a-am)>3*a+am){
    fm <- 1+gamma*(HGD(0,NIND-1,k-1,2)*(2*am/(3*(a+am))+b*am/(3*(b*(a+am)+a)))
                   +HGD(1,NIND-1,k-1,2)*b*am/(2*b*am+a)
                   +HGD(2,NIND-1,k-1,2)/3.0-c*am)
    
    fr <- 1+gamma*(HGD(0,NIND-1,k,2)/3.0
                   +HGD(1,NIND-1,k,2)*(a/(3*(a+am))+(1+b)*a/(6*(b*(a+am)+a)))
                   +HGD(2,NIND-1,k,2)*a/(2*b*am+a)-c*a)
    
  }else if(b<am*(am+3*a)/((am+a)*(am-a))){
    fm <- 1+gamma*(HGD(0,NIND-1,k-1,2)*am*(3*am+a)/(3*(am+a)*(am+a))
                   +HGD(1,NIND-1,k-1,2)*am*(3*am+5*a)/(6*(am+a)*(am+a))
                   +HGD(2,NIND-1,k-1,2)/3.0-c*am)
    
    fr <- 1+gamma*(HGD(0,NIND-1,k,2)/3.0
                   +HGD(1,NIND-1,k,2)*a*(5*am+3*a)/(6*(am+a)*(am+a))
                   +HGD(2,NIND-1,k,2)*a*(am+3*a)/(3*(am+a)*(am+a))-c*a)
    
  }else if((b>am*(am+3*a)/((am+a)*(am-a))) && (b<(3*am+a)/(2*(am-a)))){
    fm <- 1+gamma*(HGD(0,NIND-1,k-1,2)*am*(3*am+a)/(3*(am+a)*(am+a))
                   +HGD(1,NIND-1,k-1,2)*(am/(3*(am+a))+(1+b)*am/(6*(b*(am+a)+am)))
                   +HGD(2,NIND-1,k-1,2)/3.0-c*am)
    
    fr <- 1+gamma*(HGD(0,NIND-1,k,2)/3.0
                   +HGD(1,NIND-1,k,2)*a*(5*am+3*a)/(6*(am+a)*(am+a))
                   +HGD(2,NIND-1,k,2)*(2*a/(3*(am+a))+b*a/(3*(b*(am+a)+am)))-c*a)
    
  }else if(2*b*(am-a)>3*am+a){
    fm <- 1+gamma*(HGD(0,NIND-1,k-1,2)*am/(2*b*a+am)
                   +HGD(1,NIND-1,k-1,2)*(am/(3*(am+a))+(1+b)*am/(6*(b*(am+a)+am)))
                   +HGD(2,NIND-1,k-1,2)/3.0-c*am)
    
    fr <- 1+gamma*(HGD(0,NIND-1,k,2)/3.0
                   +HGD(1,NIND-1,k,2)*b*a/(2*b*a+am)
                   +HGD(2,NIND-1,k,2)*(2*a/(3*(am+a))+b*a/(3*(b*(am+a)+am)))-c*a)
    
  }else
    stop("error\n")
  return(list(fr = fr, fm = fm))
}


HGD <- function(i, n, K, m) {
  res <- choose(K,i)*choose(n-K,m-i)/choose(n,m)
  return(res)
}


Binomial <- function(s, t) {
  if(s<t) { 
    b <- 0
  }
  else{
    b <- 1
    for(i in 0:(t-1)) {  # index change might change the result
      b <- b * as.double((s-i)/(t-i+1))
    }
  }
  return(b)
}


Reproduction <- function(nm, fm, fr) {
  nmn <- 0
  for(i in 0:(NIND-1)) {
    if(runif(1)<fm*nm/(fm*nm+fr*(NIND-nm)))
      nmn <- nmn+1
  }
  return(nmn)
}




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