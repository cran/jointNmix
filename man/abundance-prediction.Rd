\name{abundance-prediction}
\alias{getranef.uniNmix}
\alias{getranef.jointNmix}
\encoding{UTF-8}
\title{Get empirical Bayes predictions of the latent abundances}
\description{Computes the conditional abundance distribution given the data and model estimates}
\usage{
getranef.uniNmix(obj, distr = FALSE)
getranef.jointNmix(obj, distr = FALSE)
}
\arguments{
  \item{obj}{fitted model object}
  \item{distr}{logical. If TRUE, returns the matrix of conditional probabilities. If FALSE, returns the expected values of the conditional distributions, i.e. the abundance estimates}
}
\details{
These functions returnt the empirical Bayes estimates of the latent abundances from univariate and joint N-mixture models, as described by Royle (2004) and Moral et al. (submitted), respectively.
}
\references{
Royle, J.A. (2004) Models for estimating population size from spatially replicated counts. Biometrics 60:108-105.

Moral, R.A., Hinde, J., \enc{Demétrio}{Demetrio}, C.G.B., Reigada, C. and Godoy, W.A.C. (submitted) Models for jointly estimating abundance of two unmarked site-associated species subject to imperfect detection.
}
\author{
Rafael A. Moral <rafael_moral@yahoo.com.br>, Clarice G. B. \enc{Demétrio}{Demetrio} and John Hinde
}

\seealso{
\code{\link{Nmix}}
\code{\link{jointNmix}}
}