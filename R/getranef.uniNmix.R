getranef.uniNmix <-
function(obj, distr=FALSE) {
  mix <- obj$mixture
  K <- obj$K
  parms <- coef(obj)
  sp1 <- obj$sp1
  np1 <- ncol(obj$Xp)
  nl1 <- ncol(obj$Xl)
  pr1 <- as.numeric(plogis(obj$Xp %*% parms[1:np1]))
  lam1 <- as.numeric(exp(obj$Xl %*% parms[(np1+1):(np1+nl1)])) # dimension R
  if(mix=="NB") theta1 <- exp(parms[(np1+nl1+1)])
  if(mix=="NeymanA") lam2 <- exp(parms[(np1+nl1+1)])
  R <- dim(sp1)[1]
  T_ <- dim(sp1)[2]
 
  ########### binomial part
  pr1 <- matrix(pr1, byrow=TRUE, ncol=T_)
  
  bin1 <- matrix(0, ncol=K+1, nrow=R)
  for(Ni1 in 0:K) {
    bin1.aux <- matrix(0, ncol=T_, nrow=R)
    for(i in 1:T_) bin1.aux[,i] <- dbinom(sp1[,i], Ni1, pr1[,i], log=TRUE)
    bin1[,Ni1+1] <- rowSums(bin1.aux)
  }

  ########### abundance part
  f1 <- matrix(0, nrow=R, ncol=K+1)
  if(mix=="P") for(i in 1:R) f1[i,] <- dpois(0:K, lam1[i], log=TRUE)
  if(mix=="NB") for(i in 1:R) f1[i,] <- dnbinom(0:K, mu=lam1[i], size=theta1, log=TRUE)
  if(mix=="NeymanA") for(i in 1:R) f1[i,] <- dneymanA(0:K, lambda1=lam1[i], lambda2=lam2, K=K, log=TRUE)
  
  num1 <- exp(bin1 + f1)
  den1 <- matrix(rowSums(num1), byrow=F, nrow=R, ncol=K+1)
  post1 <- num1/den1 * matrix(0:K, byrow=TRUE, nrow=R, ncol=K+1)
  mean1 <- round(apply(post1, 1, sum), 2)
  
  expectedN <- data.frame("Ni1"=mean1)
  if(!distr) return(expectedN) else return(list("sp1"=num1/den1))
}
