getranef.jointNmix <-
function(obj, distr=FALSE) {
  lam2transform <- obj$lam2transform
  includepsi <- obj$includepsi
  mix <- obj$mixture
  K <- obj$K
  parms <- coef(obj)
  sp1 <- obj$sp1
  sp2 <- obj$sp2
  np1 <- ncol(obj$Xp1)
  np2 <- ncol(obj$Xp2)
  nl1 <- ncol(obj$Xl1)
  nl2 <- ncol(obj$Xl2)
  npsi <- ncol(obj$Xpsi)
  pr1 <- as.numeric(plogis(obj$Xp1 %*% parms[1:np1]))
  pr2 <- as.numeric(plogis(obj$Xp2 %*% parms[(np1+1):(np1+np2)]))
  lam1 <- as.numeric(exp(obj$Xl1 %*% parms[(np1+np2+1):(np1+np2+nl1)]))
  lam2 <- as.numeric(exp(obj$Xl2 %*% parms[(np1+np2+nl1+1):(np1+np2+nl1+nl2)]))
  ########## lam2 transform
  if(lam2transform) lam2 <- lam2/(1+lam2)
  ########## includepsi
  if(includepsi) psi <- as.numeric(exp(obj$Xpsi %*% parms[(length(parms)-npsi+1):(length(parms))]))
  ########## negbin theta (phi)
  if(mix[1]=="NB") theta1 <- exp(parms[(np1+np2+nl1+nl2+1)])
  if(mix[2]=="NB"&mix[1]=="P") theta2 <- exp(parms[(np1+np2+nl1+nl2+1)])
  if(mix[2]=="NB"&mix[1]=="NB") theta2 <- exp(parms[(np1+np2+nl1+nl2+2)])
  
  R <- dim(sp1)[1]
  T_ <- dim(sp1)[2]
  
  ###########
  pr1 <- matrix(pr1, byrow=TRUE, ncol=T_)
  pr2 <- matrix(pr2, byrow=TRUE, ncol=T_)
  
  bin1 <- bin2 <- matrix(0, ncol=K+1, nrow=R)
  for(Ni1 in 0:K) {
    bin1.aux <- matrix(0, ncol=T_, nrow=R)
    for(i in 1:T_) bin1.aux[,i] <- dbinom(sp1[,i], Ni1, pr1[,i], log=TRUE)
    bin1[,Ni1+1] <- rowSums(bin1.aux)
  }
  for(Ni2 in 0:K) {
    bin2.aux <- matrix(0, ncol=T_, nrow=R)
    for(i in 1:T_) bin2.aux[,i] <- dbinom(sp2[,i], Ni2, pr2[,i], log=TRUE)
    bin2[,Ni2+1] <- rowSums(bin2.aux)
  }
  ###########
  
  ###########
  f1 <- matrix(0, nrow=R, ncol=K+1)
  if(mix[1]=="P") for(i in 1:R) f1[i,] <- dpois(0:K, lam1[i], log=TRUE)
  if(mix[1]=="NB") for(i in 1:R) f1[i,] <- dnbinom(0:K, mu=lam1[i], size=theta1, log=TRUE)
  ###########
  
  num1 <- exp(bin1 + f1)
  den1 <- matrix(rowSums(num1), byrow=F, nrow=R, ncol=K+1)
  post1 <- num1/den1 * matrix(0:K, byrow=TRUE, nrow=R, ncol=K+1)
  mean1 <- round(apply(post1, 1, sum), 2)
  
  ###########
  f2 <- matrix(0, nrow=R, ncol=K+1)
  for(Ni2 in 0:K) {
    f2.aux2 <- matrix(0, nrow=R, ncol=K+1)
    for(Ni1 in 0:K) {
      f2.aux <- NULL
      if(mix[1]=="P"&mix[2]=="P") {
        if(!includepsi) {
          for(i in 1:R) f2.aux[i] <- dpois(Ni1, lam1[i]) * dpois(Ni2, lam2[i] * Ni1)
        } else {
          for(i in 1:R) f2.aux[i] <- dpois(Ni1, lam1[i]) * dpois(Ni2, psi[i] + lam2[i] * Ni1)
        }
      }
      if(mix[1]=="P"&mix[2]=="NB") {
        if(!includepsi) {
          for(i in 1:R) for(i in 1:R) f2.aux[i] <- dpois(Ni1, lam1[i]) * dnbinom(Ni2, mu=lam2[i] * Ni1, size=theta2)
        } else {
          for(i in 1:R) for(i in 1:R) f2.aux[i] <- dpois(Ni1, lam1[i]) * dnbinom(Ni2, mu=psi[i] + lam2[i] * Ni1, size=theta2)
        }
      }
      if(mix[1]=="NB"&mix[2]=="P") {
        if(!includepsi) {
          for(i in 1:R) f2.aux[i] <- dnbinom(Ni1, mu=lam1[i], size=theta1) * dpois(Ni2, lam2[i] * Ni1)
        } else {
          for(i in 1:R) f2.aux[i] <- dnbinom(Ni1, mu=lam1[i], size=theta1) * dpois(Ni2, psi[i] + lam2[i] * Ni1)
        }
      }
      if(mix[1]=="NB"&mix[2]=="NB") {
        if(!includepsi) {
          for(i in 1:R) f2.aux[i] <- dnbinom(Ni1, mu=lam1[i], size=theta1) * dnbinom(Ni2, mu=lam2[i] * Ni1, size=theta2)
        } else {
          for(i in 1:R) f2.aux[i] <- dnbinom(Ni1, mu=lam1[i], size=theta1) * dnbinom(Ni2, mu=psi[i] + lam2[i] * Ni1, size=theta2)
        }
      }
      f2.aux2[,Ni1+1] <- f2.aux
    }
    f2[,Ni2+1] <- rowSums(f2.aux2)
  }
  ###########
  
  num2 <- exp(bin2) * f2
  den2 <- matrix(rowSums(num2), byrow=FALSE, nrow=R, ncol=K+1)
  post2 <- num2/den2 * matrix(0:K, byrow=TRUE, nrow=R, ncol=K+1)
  mean2 <- round(apply(post2, 1, sum), 2)
  expectedN <- data.frame("Ni1"=mean1, "Ni2"=mean2)
  if(!distr) return(expectedN) else return(list("sp1"=num1/den1, "sp2"=num2/den2))
}
