jointNmix <-
function(sp1, sp2, start, method="BFGS", K, mixture=c("P","P"), 
                      Xp1, Xp2, Xl1, Xl2, Xpsi, includepsi=TRUE) {
  # Xl1, Xl2 and Xpsi are design matrices (R x p)
  # Xp1 and Xp2 are design matrices (RT x p)
  if(missing(mixture)) mixture <- c("P","P")
  if(missing(Xpsi)) Xpsi <- Xl2 
  if(missing(start)) {
    if(mixture[1]=="P"&mixture[2]=="P") start <- rep(0, ncol(Xp1)+ncol(Xp2)+ncol(Xl1)+ncol(Xl2))
    if(mixture[1]=="P"&mixture[2]=="NB") start <- rep(0, ncol(Xp1)+ncol(Xp2)+ncol(Xl1)+ncol(Xl2)+1)
    if(mixture[1]=="NB"&mixture[2]=="P") start <- rep(0, ncol(Xp1)+ncol(Xp2)+ncol(Xl1)+ncol(Xl2)+1)
    if(mixture[1]=="NB"&mixture[2]=="NB") start <- rep(0, ncol(Xp1)+ncol(Xp2)+ncol(Xl1)+ncol(Xl2)+2)
  }
  if(includepsi) start <- c(start, rep(0, ncol(Xpsi))) else psi.e <- 0
  if(missing(K)) {
    K <- max(sp1, sp2) + 100
    cat("K set as", K, "\n")
  }
  np1 <- ncol(Xp1)
  np2 <- ncol(Xp2)
  nl1 <- ncol(Xl1)
  nl2 <- ncol(Xl2)
  npsi <- ncol(Xpsi)
  R <- dim(sp1)[1]
  T_ <- dim(sp1)[2]
  
  #### likelihood function to be maximized
  lik <- function(parms, obj1, obj2, K, Xp1, Xp2, Xl1, Xl2, Xpsi, mixture, includepsi) {
    pr1 <- as.numeric(plogis(Xp1 %*% parms[1:np1])) # dimension RT
    pr2 <- as.numeric(plogis(Xp2 %*% parms[(np1+1):(np1+np2)])) # dimension RT
    lam1 <- as.numeric(exp(Xl1 %*% parms[(np1+np2+1):(np1+np2+nl1)])) # dimension R
    lam2 <- as.numeric(exp(Xl2 %*% parms[(np1+np2+nl1+1):(np1+np2+nl1+nl2)])) # dimension R
    ########## includepsi
    if(includepsi) psi <- as.numeric(exp(Xpsi %*% parms[(length(parms)-npsi+1):length(parms)])) # dimension R
    ########## negbin theta (phi)
    if(mixture[1]=="NB") theta1 <- exp(parms[(np1+np2+nl1+nl2+1)])
    if(mixture[2]=="NB"&mixture[1]=="P") theta2 <- exp(parms[(np1+np2+nl1+nl2+1)])
    if(mixture[2]=="NB"&mixture[1]=="NB") theta2 <- exp(parms[(np1+np2+nl1+nl2+2)])
    
    ########### binomial likelihood
    pr1 <- matrix(pr1, byrow=TRUE, ncol=T_)
    pr2 <- matrix(pr2, byrow=TRUE, ncol=T_)
    
    bin1 <- bin2 <- matrix(0, ncol=K+1, nrow=R)
    for(Ni1 in 0:K) {
      bin1.aux <- matrix(0, ncol=T_, nrow=R)
      for(i in 1:T_) bin1.aux[,i] <- dbinom(obj1[,i], Ni1, pr1[,i], log=TRUE)
      bin1[,Ni1+1] <- exp(rowSums(bin1.aux))
    }
    for(Ni2 in 0:K) {
      bin2.aux <- matrix(0, ncol=T_, nrow=R)
      for(i in 1:T_) bin2.aux[,i] <- dbinom(obj2[,i], Ni2, pr2[,i], log=TRUE)
      bin2[,Ni2+1] <- exp(rowSums(bin2.aux))
    }
    
    ########### abundance likelihood for first species
    f1 <- matrix(0, nrow=R, ncol=K+1)
    if(mixture[1]=="P") for(i in 1:R) f1[i,] <- dpois(0:K, lam1[i])
    if(mixture[1]=="NB") for(i in 1:R) f1[i,] <- dnbinom(0:K, mu=lam1[i], size=theta1)
    
    part1 <- rowSums(bin1 * f1)
    
    ########### abundance likelihood for second species
    f2 <- list()
    for(Ni1 in 0:K) {
      f2.aux <- matrix(0, nrow=R, ncol=K+1)
      if(mixture[2]=="P") {
        if(!includepsi) {
          for(i in 1:R) f2.aux[i,] <- dpois(0:K, lam2[i] * Ni1)
        } else {
          for(i in 1:R) f2.aux[i,] <- dpois(0:K, psi[i] + lam2[i] * Ni1)
        }
      }
      if(mixture[2]=="NB") {
        if(!includepsi) {
          for(i in 1:R) f2.aux[i,] <- dnbinom(0:K, mu=lam2[i] * Ni1, size=theta2)
        } else {
          for(i in 1:R) f2.aux[i,] <- dnbinom(0:K, mu=psi[i] + lam2[i] * Ni1, size=theta2)
        }
      }
      f2[[Ni1+1]] <- f2.aux
    }
    
    ########### combining
    sum1 <- lapply(f2, function(x) rowSums(x * bin2))
    sum1 <- as.matrix(data.frame(sum1))
    colnames(sum1) <- NULL
    part2 <- rowSums(sum1 * f1)
    
    ########## total likelihood
    totalsum <- part1 * part2
    loglik <- log(totalsum)
    loglik[loglik==-Inf] <- -9999
    return(-sum(loglik))
  }
  fit.obj <- optim(start, lik, obj1=sp1, obj2=sp2, hessian=TRUE, method=method, K=K, Xp1=Xp1, Xp2=Xp2, Xl1=Xl1, Xl2=Xl2, Xpsi=Xpsi, mixture=mixture, includepsi=includepsi)
  fit.pars <- fit.obj$par
  names(fit.pars) <- c(paste("p1:", colnames(Xp1)), paste("p2:", colnames(Xp2)), paste("lambda1:", colnames(Xl1)), paste("lambda2:", colnames(Xl2)))
  if(mixture[1]=="NB"&mixture[2]=="P") names(fit.pars)[np1+np2+nl1+nl2+1] <- "theta1"
  if(mixture[1]=="P"&mixture[2]=="NB") names(fit.pars)[np1+np2+nl1+nl2+1] <- "theta2"
  if(mixture[1]=="NB"&mixture[2]=="NB") names(fit.pars)[(np1+np2+nl1+nl2+1):(np1+np2+nl1+nl2+2)] <- c("theta1","theta2")
  if(includepsi) names(fit.pars)[(length(fit.pars)-npsi+1):(length(fit.pars))] <- paste("psi:", colnames(Xl2))
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
  p1.e <- plogis(Xp1 %*% fit.pars[1:np1])
  p2.e <- plogis(Xp2 %*% fit.pars[(np1+1):(np1+np2)])
  p1.e <- matrix(p1.e, byrow=TRUE, ncol=T_)
  p2.e <- matrix(p2.e, byrow=TRUE, ncol=T_)
  l1.e <- exp(Xl1 %*% fit.pars[(np1+np2+1):(np1+np2+nl1)])
  l2.e <- exp(Xl2 %*% fit.pars[(np1+np2+nl1+1):(np1+np2+nl1+nl2)])
  ########### includepsi
  if(includepsi) psi.e <- exp(Xpsi %*% fit.pars[(length(fit.pars)-npsi+1):length(fit.pars)])
  if(mixture[1]=="NB") t1.e <- exp(fit.pars[(np1+np2+nl1+nl2+1)])
  if(mixture[2]=="NB"&mixture[1]=="P") t2.e <- exp(fit.pars[(np1+np2+nl1+nl2+1)])
  if(mixture[2]=="NB"&mixture[1]=="NB") t2.e <- exp(fit.pars[(np1+np2+nl1+nl2+2)])
  
  #### fitted values
  fv1 <- matrix(l1.e, byrow=F, ncol=T_, nrow=R) * p1.e
  if(!includepsi) {
    fv2 <- matrix(l1.e, byrow=F, ncol=T_, nrow=R) * matrix(l2.e, byrow=F, ncol=T_, nrow=R) * p2.e
  } else {
    fv2 <- (matrix(l1.e, byrow=F, ncol=T_, nrow=R) * matrix(l2.e, byrow=F, ncol=T_, nrow=R) + matrix(psi.e, byrow=F, ncol=T_, nrow=R))* p2.e
  }
  
  #### weights for residuals
  if(mixture[1]=="P") weights1 <- fv1 else weights1 <- fv1 * (1 + 1/t1.e)
  if(mixture[2]=="P") weights2 <- fv2 else weights2 <- fv2 * (1 + 1/t2.e)
  
  #### covariance and correlation - observed Y
  covar <- corr <- matrix(0, ncol=T_, nrow=R)
  if(mixture[1]=="P"&mixture[2]=="P") for(i in 1:T_) covar[,i] <- l1.e*l2.e*p1.e[,i]*p2.e[,i]
  if(mixture[1]=="P"&mixture[2]=="P") for(i in 1:T_) corr[,i] <- l2.e*sqrt(l1.e*p1.e[,i]*p2.e[,i]/(psi.e+l1.e*l2.e*(1+l2.e*p2.e[,i])))
  if(mixture[1]=="P"&mixture[2]=="NB") for(i in 1:T_) covar[,i] <- l1.e*l2.e*p1.e[,i]*p2.e[,i]
  if(mixture[1]=="P"&mixture[2]=="NB") for(i in 1:T_) corr[,i] <- l2.e*sqrt(p1.e[,i]*p2.e[,i]*l1.e/(psi.e+l1.e*l2.e*(1+l2.e*p2.e[,i])+p2.e[,i]/t2.e*(l1.e*l2.e^2+(psi.e+l1.e*l2.e)^2)))
  if(mixture[1]=="NB"&mixture[2]=="P") for(i in 1:T_) covar[,i] <- l1.e*l2.e*p1.e[,i]*p2.e[,i]*(l1.e^2/t1.e)
  if(mixture[1]=="NB"&mixture[2]=="P") for(i in 1:T_) corr[,i] <- l2.e*(1+l1.e/t1.e)*sqrt(l1.e*p1.e[,i]*p2.e[,i]/((1+p1.e[,i]*l1.e/t1.e)*(psi.e+l1.e*l2.e*(1+l2.e*p2.e[,i]*(1+l1.e/t1.e)))))
  if(mixture[1]=="NB"&mixture[2]=="NB") for(i in 1:T_) covar[,i] <- l1.e*l2.e*p1.e[,i]*p2.e[,i]*(l1.e^2/t1.e)
  if(mixture[1]=="NB"&mixture[2]=="NB") for(i in 1:T_) corr[,i] <- l2.e*(1+l1.e/t1.e)*sqrt(l1.e*p1.e[,i]*p2.e[,i]/((1+p1.e[,i]*l1.e/t1.e)*(psi.e+l1.e*l2.e*(1+l2.e*p2.e[,i]*(1+l1.e/t1.e))+l1.e*l2.e*p2.e[,i]/(t1.e*t2.e))+p2.e[,i]/t2.e*(l1.e*l2.e^2+(psi.e+l1.e*l2.e)^2)))
    
  Nmix.result <- list("coef"=fit.pars, "se"=fit.se, "z"=fit.z, "p"=fit.p, 
                      "logLik"=llik, "AIC"=AIC, "K"=K, "sp1"=sp1, "sp2"=sp2, 
                      "Xp1"=Xp1, "Xp2"=Xp2, "Xl1"=Xl1, "Xl2"=Xl2, "Xpsi"=Xpsi,
                      "mixture"=mixture, 
                      "Cov"=covar, "Corr"=corr, "fv1"=fv1, "fv2"=fv2, "weights1"=weights1, "weights2"=weights2,
                      "lam2transform"=FALSE, "includepsi"=includepsi)
  class(Nmix.result) <- c("jointNmix","Nmix")
  cat("Joint N-mixture model fit by maximum likelihood", "\n")
  return(Nmix.result)
}
