\name{Nmix}
\alias{Nmix}
\encoding{UTF-8}
\title{
Univariate N-mixture models
}
\description{
Fits univariate N-mixture models
}
\usage{
Nmix(sp1, start, method = "BFGS", K, mixture, Xp, Xl)
}
\arguments{
  \item{sp1}{observation matrix for the species}
  \item{start}{initial values for the optimization process}
  \item{method}{optimization method passed to \code{\link{optim}}. Defaults to \code{"BFGS"}}
  \item{K}{truncation number of the infinite summations in the log-likelihood. Defaults to \code{max(sp1, sp2) + 100}}
  \item{mixture}{latent abundance distribution specification. \code{"P"} for Poisson, \code{"NB"} for negative binomial and \code{"NeymanA"} for Neyman-A distributions are available. Defaults to \code{"P"}}
  \item{Xp}{model matrix for detection probabilities}
  \item{Xl}{model matrix for abundance parameter}
}
\details{
The function fits Royle's (2004) N-mixture model to data on species abundance collected at R sites over T time occasions. The model for observation on site i at time t can be specified as \deqn{Y_{it}|N_i ~ Bin(N_i,p_{it})}\deqn{N_i ~ a count distribution with mean \lambda_i.} Here, users may define a Poisson, negative binomial or Neyman-A distributions for the latent abundances N_i.
}
\value{
An object of class \code{uniNmix} and \code{Nmix}, for which many methods are available (see \code{methods(class = "uniNmix")} and \code{methods(class = "Nmix")})
}
\references{
Royle, J.A. (2004) Models for estimating population size from spatially replicated counts. Biometrics 60:108-105.
}
\author{
Rafael A. Moral <rafael_moral@yahoo.com.br>, Clarice G. B. \enc{Demétrio}{Demetrio} and John Hinde
}\seealso{
\code{\link{jointNmix}}
}
\examples{
## simulating data with negative binomial latent abundances
R <- 10 # sites
T <- 10 # time occasions
lambda <- 5 # abundance parameter
p <- .3 # probability of detection
phi <- 1 # dispersion parameter
set.seed(1234); Ni <- rnbinom(R, mu=lambda, size=phi) # latent abundances
y <- matrix(0, ncol=T, nrow=R)
set.seed(1234); for(i in 1:R) y[,i] <- rbinom(T, Ni, p) # observed abundances

## fitting the Poisson N-mixture model
fitp <- Nmix(y, Xp=cbind(rep(1, R*T)), Xl=cbind(rep(1, R)), mixture="P", K=25)

## fitting the negative binomial N-mixture model
fitnb <- Nmix(y, Xp=cbind(rep(1, R*T)), Xl=cbind(rep(1, R)), mixture="NB", K=25)

## fitting the Neyman-A N-mixture model
fitna <- Nmix(y, Xp=cbind(rep(1, R*T)), Xl=cbind(rep(1, R)), mixture="NeymanA", K=25)

## likelihood-ratio test between Poisson and negbin models
anova(fitp, fitnb)

## comparing using AIC
lapply(list(fitp, fitnb, fitna), AIC)

## conditional posterior probability functions for abundances
plot(fitnb, posterior = TRUE)

## estimated abundances vs. true abundances
data.frame(getranef.uniNmix(fitnb), Ni)
}