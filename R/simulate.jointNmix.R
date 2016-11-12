simulate.jointNmix <-
function(object, ...) {
  lam2transform <- object$lam2transform
  includepsi <- object$includepsi
  parms <- coef(object)
  np1 <- ncol(object$Xp1)
  np2 <- ncol(object$Xp2)
  nl1 <- ncol(object$Xl1)
  nl2 <- ncol(object$Xl2)
  npsi <- ncol(object$Xpsi)
  pr1 <- as.numeric(plogis(object$Xp1 %*% parms[1:np1])) # dimension RT
  pr2 <- as.numeric(plogis(object$Xp2 %*% parms[(np1+1):(np1+np2)])) # dimension RT
  lam1 <- as.numeric(exp(object$Xl1 %*% parms[(np1+np2+1):(np1+np2+nl1)])) # dimension R
  lam2 <- as.numeric(exp(object$Xl2 %*% parms[(np1+np2+nl1+1):(np1+np2+nl1+nl2)])) # dimension R
  ########## lam2 transform
  if(lam2transform) lam2 <- lam2/(1+lam2)
  ########## includepsi
  if(includepsi) psi <- as.numeric(exp(object$Xpsi %*% parms[(length(parms)-npsi+1):(length(parms))]))
  theta1 <- exp(parms[(np1+np2+nl1+nl2+1)])
  theta2 <- exp(parms[length(parms)])

  R <- dim(object$sp1)[1]
  T_ <- dim(object$sp1)[2]
  
  pr1 <- matrix(pr1, byrow=TRUE, ncol=T_)
  pr2 <- matrix(pr2, byrow=TRUE, ncol=T_)
  
  ## simulating the latent abundances
  ### true abundance of sp. 1
  if(object$mixture[1]=="P") Ni1 <- rpois(R, lam1)
  if(object$mixture[1]=="NB") Ni1 <- rnbinom(R, mu=lam1, size=theta1) 
  ### true abundance of sp. 2
  if(object$mixture[2]=="P") {
    if(!includepsi) {
      Ni2 <- rpois(R, lam2 * Ni1)
    } else {
      Ni2 <- rpois(R, psi + lam2 * Ni1)
    }
  }
  if(object$mixture[2]=="NB") {
    if(!includepsi) {
      Ni2 <- rnbinom(R, mu=lam2 * Ni1, size=theta2)
    } else {
      Ni2 <- rnbinom(R, mu=psi + lam2 * Ni1, size=theta2)
    }
  }
  
  ## simulating the observations
  sdata1 <- matrix(0, ncol=T_, nrow=R)
  sdata2 <- matrix(0, ncol=T_, nrow=R)
  for(i in 1:T_) sdata1[,i] <- rbinom(R, Ni1, pr1[,i])
  for(i in 1:T_) sdata2[,i] <- rbinom(R, Ni2, pr2[,i])
  
  return(list("sp1"=sdata1, "sp2"=sdata2))
}
