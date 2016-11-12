simulate.uniNmix <-
function(object, ...) {
  parms <- coef(object)
  np1 <- ncol(object$Xp)
  nl1 <- ncol(object$Xl)
  pr1 <- as.numeric(plogis(object$Xp %*% parms[1:np1])) # dimension RT
  lam1 <- as.numeric(exp(object$Xl %*% parms[(np1+1):(np1+nl1)])) # dimension R
  if(object$mixture=="NB") theta1 <- exp(parms[(np1+nl1+1)])
  if(object$mixture=="NeymanA") lam2 <- exp(parms[(np1+nl1+1)])
  R <- dim(object$sp1)[1]
  T_ <- dim(object$sp1)[2]
  
  pr1 <- matrix(pr1, byrow=TRUE, ncol=T_)
  
  ## simulating the latent abundances
  ### true abundance
  if(object$mixture=="P") Ni1 <- rpois(R, lam1)
  if(object$mixture=="NB") Ni1 <- rnbinom(R, mu=lam1, size=theta1) 
  if(object$mixture=="NeymanA") {
    NN <- rpois(R, lam1)
    Ni1 <- rpois(R, lam2 * NN)
  }
  
  ## simulating the observations
  sdata1 <- matrix(0, ncol=T_, nrow=R)
  for(i in 1:T_) sdata1[,i] <- rbinom(R, Ni1, pr1[,i])
  
  return(list("sp1"=sdata1))
}
