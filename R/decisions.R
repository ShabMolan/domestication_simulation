source("R/functions.R")

c <- 0.05
ar <- 5

aS <- Mutation(ar, c)
aT <- Mutation(ar, c)
aH <- Mutation(ar, c)

b <- 5

alpha <- 1000
beta <- 1000

wSS <- (b * aS) / (b * (aS + aH) + aT)
wTS <- (aT) / (b * (aS + aH) + aT)
wHS <- (b * aH) / (b * (aS + aH) + aT)

wST <- (aS) / (b * (aT + aH) + aS)
wTT <- (b * aT) / (b * (aT + aH) + aS)
wHT <- (b * aH) / (b * (aT + aH) + aS)

wSN <- (aS ^ 2) / ((aS + aT) * (aS + aH))
wTN <- (aT ^ 2) / ((aS + aT) * (aS + aH))
wHN <- ((aH * aS) / ((aS + aT) * (aS + aH))) + 
       ((aH * aT) / ((aS + aT) * (aS + aH)))

g2 <- (1 / (exp(alpha * wHS) + exp(alpha * wHT) + exp(alpha * wHN))) * 
      c(exp(alpha * wHS), exp(alpha * wHT), exp(alpha * wHN))

Helper_decision <- Helper_decision(g2)

wS <- Helper_decision * c(wSS, wST, wSN)

g1 <- (1 / (exp(beta * wS) + exp(beta * wS))) * 
      c(exp(beta * wS), exp(beta * wS))



  





