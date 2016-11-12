Nmix <-
function(sp1, start, method="BFGS", K, mixture, Xp, Xl) {
  if(missing(mixture)) mixture <- "P"
  if(missing(start)) {
    if(mixture=="P") start <- rep(0, ncol(Xp)+ncol(Xl)) else if(mixture=="NB"|mixture=="NeymanA") start <- rep(0, ncol(Xp)+ncol(Xl)+1) else stop("mixture must be 'P' for Poisson, 'NB' for negative binomial or 'NA' for Neyman-A")
  }
  if(missing(K)) {
    K <- max(sp1) + 100
    cat("K set as", K, "\n")
  }
  np1 <- ncol(Xp)
  nl1 <- ncol(Xl)
  R <- dim(sp1)[1]
  T_ <- dim(sp1)[2]
  
  #### likelihood function to be maximized
  lik <- function(parms, obj1, K, Xp, Xl, mixture) {
    pr1 <- as.numeric(plogis(Xp %*% parms[1:np1])) # dimension RT
    lam1 <- as.numeric(exp(Xl %*% parms[(np1+1):(np1+nl1)])) # dimension R
    if(mixture=="NB") theta1 <- exp(parms[(np1+nl1+1)])
    if(mixture=="NeymanA") lam2 <- exp(parms[(np1+nl1+1)])
    
    ########### binomial likelihood
    pr1 <- matrix(pr1, byrow=TRUE, ncol=T_)

    bin1 <- matrix(0, ncol=K+1, nrow=R)
    for(Ni1 in 0:K) {
      bin1.aux <- matrix(0, ncol=T_, nrow=R)
      for(i in 1:T_) bin1.aux[,i] <- dbinom(obj1[,i], Ni1, pr1[,i], log=TRUE)
      bin1[,Ni1+1] <- exp(rowSums(bin1.aux))
    }
    
    ########### abundance likelihood
    f1 <- matrix(0, nrow=R, ncol=K+1)
    if(mixture=="P") for(i in 1:R) f1[i,] <- dpois(0:K, lam1[i])
    if(mixture=="NB") for(i in 1:R) f1[i,] <- dnbinom(0:K, mu=lam1[i], size=theta1)
    if(mixture=="NeymanA") for(i in 1:R) f1[i,] <- dneymanA(0:K, lambda1=lam1[i], lambda2=lam2, K=K)
    
    ########## total likelihood
    totalsum <- rowSums(bin1 * f1)
    loglik <- log(totalsum)
    loglik[loglik==-Inf] <- -9999
    return(-sum(loglik))
  }
  fit.obj <- optim(start, lik, obj1=sp1, hessian=TRUE, method=method, K=K, Xp=Xp, Xl=Xl, mixture=mixture)
  fit.pars <- fit.obj$par
  names(fit.pars) <- c(paste("p:", colnames(Xp)), paste("lambda:", colnames(Xl)))
  if(mixture=="NB") names(fit.pars)[np1+nl1+1] <- "theta"
  if(mixture=="NeymanA") names(fit.pars)[np1+nl1+1] <- "lambda2"
  fit.se <- try(diag(solve(fit.obj$hessian)), silent=TRUE)
  if(class(fit.se)=="try-error") {
    warning("The Hessian matrix was singular; se's not computed")
    fit.se <- fit.z <- fit.p <- NA
  } else {
    fit.z <- fit.pars/fit.se
    fit.p <- pnorm(abs(fit.z), lower.tail=F)*2
  }
  llik <- fit.obj$value
  AIC <- 2*(llik + length(fit.pars))
  p1.e <- plogis(Xp %*% fit.pars[1:np1])
  p1.e <- matrix(p1.e, byrow=TRUE, ncol=T_)
  l1.e <- exp(Xl %*% fit.pars[(np1+1):(np1+nl1)])
  if(mixture=="NB") t1.e <- exp(fit.pars[(np1+nl1+1)])
  if(mixture=="NeymanA") l2.e <- exp(fit.pars[(np1+nl1+1)])
  
  #### fitted values
  if(mixture=="P"|mixture=="NB") fv1 <- matrix(l1.e, byrow=F, ncol=T_, nrow=R) * p1.e
  if(mixture=="NeymanA") fv1 <- matrix(l1.e, byrow=F, ncol=T_, nrow=R) * p1.e * l2.e
  
  #### weights for residuals
  if(mixture=="P") weights1 <- fv1 else if(mixture=="NB") weights1 <- fv1 * (1 + 1/t1.e) else if(mixture=="NeymanA") weights1 <- fv1 * (1 + l2.e)
  
  Nmix.result <- list("coef"=fit.pars, "se"=fit.se, "z"=fit.z, "p"=fit.p, "logLik"=llik, "AIC"=AIC, "K"=K, "sp1"=sp1, "Xp"=Xp, "Xl"=Xl, "mixture"=mixture, "fv1"=fv1, "weights1"=weights1)
  class(Nmix.result) <- c("uniNmix", "Nmix")
  cat("Univariate N-mixture model fit by maximum likelihood", "\n")
  return(Nmix.result)
}
