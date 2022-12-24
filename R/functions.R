VERS <- 2.00
SEED <- 13
NIND <- 1000
TMAX <- 1100  # 110000
BURN <- 100   # 10000
NRUN <- 10  # 100


#' Title
#' 
#' 
#'
#' @param ar 
#' @param c 
#'
#' @return
Mutation <- function(ar, c) { # why c and what is ar
  return(2*1.0/(4*c) * runif(1))
}


#' Title
#'
#' @param b 
#' @param c 
#' @param gamma 
#' @param am 
#' @param a 
#' @param k 
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param i 
#' @param n 
#' @param K 
#' @param m 
#'
#' @return
#' @export
#'
#' @examples
HGD <- function(i, n, K, m) {
  res <- choose(K,i)*choose(n-K,m-i)/choose(n,m)
  return(res)
}




#' Title
#'
#' @param nm 
#' @param fm 
#' @param fr 
#'
#' @return
#' @export
#'
#' @examples
Reproduction <- function(nm, fm, fr) {
  nmn <- 0
  for(i in 0:(NIND-1)) {
    if(runif(1)<fm*nm/(fm*nm+fr*(NIND-nm)))
      nmn <- nmn+1
  }
  return(nmn)
}
