VERS <- 2.00
SEED <- 123
NIND <- 1000
TMAX <- 1100  # 110000
BURN <- 100   # 10000
NRUN <- 10  # 100



#' Mutation
#' @param ar resident ability
#' @param c cost of ability
#'
#' @return am (mutant ability)
Mutation <- function(ar, c) {
  return(ar/(2*c) * runif(1))
}



#' Hyper Geometric Distribution (HGD)
#' @param i sample size of mutants
#' @param n population size
#' @param k number of mutants
#' @param m sample size of residents
#'
#' @return probability of drawing mutants
HGD <- function(i, n, k, m) {
  res <- (choose(k , i) * choose(n-k , m-i)) / choose(n , m)
  return(res)
}



#' Fitness
#'
#' @param b synergy
#' @param c cost of ability
#' @param gamma selection strength
#' @param am mutant ability
#' @param a ability
#' @param k number of mutants
#'
#' @return fitness
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



#' Reproduction
#' @param nm number of mutants
#' @param fm mutant fitness
#' @param fr resident fitness
#'
#' @return number of mutants in next generation
Reproduction <- function(nm, fm, fr) {
  nmn <- 0
  for(i in 0:(NIND-1)) {
    if(runif(1)<fm*nm/(fm*nm+fr*(NIND-nm)))
      nmn <- nmn+1
  }
  return(nmn)
}



Helper_decision <- function(g2) {
  if (g2[1] > g2[2] && g2[1] > g2[3]) {g2 <- c(1, 0, 0)}
  else if (g2[2] > g2[1] && g2[2] > g2[3]) {g2 <- c(0, 1, 0)}
  else {g2 <- c(0, 0, 1)}
  return(g2)
}

