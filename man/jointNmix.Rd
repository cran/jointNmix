\name{jointNmix}
\alias{jointNmix}
\encoding{UTF-8}
\title{
Joint N-mixture models 
}
\description{
Fits joint N-mixture models for site-associated species
}
\usage{
jointNmix(sp1, sp2, start, method = "BFGS", K, mixture = c("P", "P"), 
          Xp1, Xp2, Xl1, Xl2, Xpsi, includepsi = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp1}{observation matrix for species 1}
  \item{sp2}{observation matrix for species 2}
  \item{start}{initial values for the optimization process}
  \item{method}{optimization method passed to \code{\link{optim}}. Defaults to \code{"BFGS"}}
  \item{K}{truncation number of the infinite summations in the log-likelihood. Defaults to \code{max(sp1, sp2) + 100}}
  \item{mixture}{two-character vector for latent abundance distributions. \code{"P"} for Poisson and \code{"NB"} for negative binomial. Defaults to \code{c("P","P")}}
  \item{Xp1}{model matrix for detection probabilities of species 1}
  \item{Xp2}{model matrix for detection probabilities of species 2}
  \item{Xl1}{model matrix for abundance of species 1}
  \item{Xl2}{model matrix for linking parameter of species 2}
  \item{Xpsi}{model matrix for abundance of species 2}
  \item{includepsi}{logical. If FALSE, psi is not estimated and set to zero}
}
\details{
The function fits a bivariate extension to Royle's (2004) N-mixture model to data on the abundance of two species collected at R sites over T time occasions. The model for observation on site i at time t for species 1 can be specified as \deqn{Y_{1it}|N_{1i} ~ Bin(N_{1i},p_{1it})}\deqn{N_{1i} ~ a count distribution with mean \lambda_{1i}.} The model for species 2 is \deqn{Y_{2it}|N_{1i},N_{2i} ~ Bin(N_{2i},p_{2it})}\deqn{N_{2i}|N_{1i} ~ a count distribution with mean \psi+\lambda_{2i}N_{1i}.} Here, users may define a Poisson or negative binomial distribution for the latent abundances N_{1i} and N_{2i}.
}
\value{
An object of class \code{jointNmix} and \code{Nmix}, for which many methods are available (see \code{methods(class = "jointNmix")} and \code{methods(class = "Nmix")})
}
\references{Moral, R.A., Hinde, J., \enc{Demétrio}{Demetrio}, C.G.B., Reigada, C. and Godoy, W.A.C. (submitted) Models for jointly estimating abundance of two unmarked site-associated species subject to imperfect detection.}
\author{
Rafael A. Moral <rafael_moral@yahoo.com.br>, Clarice G. B. \enc{Demétrio}{Demetrio} and John Hinde
}
\seealso{
\code{\link{Nmix}}
}
\examples{
## simulating data with poisson latent abundances
R <- 10 # sites
T <- 10 # time occasions
lambda1 <- 5
psi <- 3
p1 <- .3
p2 <- .6
lambda2 <- .5
set.seed(1234); N1 <- rpois(R, lambda1)
set.seed(1234); N2 <- rpois(R, psi + lambda2*N1)
y1 <- y2 <- matrix(0, ncol=T, nrow=R)
set.seed(1234); for(i in 1:R) y1[,i] <- rbinom(T, N1, p1)
set.seed(1234); for(i in 1:R) y2[,i] <- rbinom(T, N2, p2)

Xp <- cbind(rep(1, R*T))
Xl <- cbind(rep(1, R))

\dontrun{
## fitting the Poisson-Poisson joint N-mixture model
fitpp <- jointNmix(y1, y2, Xp1=Xp, Xp2=Xp, Xl1=Xl, Xl2=Xl, mixture=c("P","P"), K=30)

## fitting the negbin-Poisson joint N-mixture model
fitnbp <- jointNmix(y1, y2, Xp1=Xp, Xp2=Xp, Xl1=Xl, Xl2=Xl, mixture=c("NB","P"), K=30)

## likelihood-ratio test between P-P and NB-P models
anova(fitpp, fitnbp)

## comparing using AIC
lapply(list(fitpp, fitnbp), AIC)

## conditional posterior probability functions for abundances
plot(fitpp, posterior = TRUE)

## estimated abundances vs. true abundances
data.frame(getranef.jointNmix(fitpp), N1, N2)
  }
}