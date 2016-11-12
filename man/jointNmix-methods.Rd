\name{jointNmix-methods}
\alias{fitted.jointNmix}
\alias{print.jointNmix}
\alias{residuals.jointNmix}
\alias{simulate.jointNmix}
\encoding{UTF-8}
\title{Methods for joint Nmix objects}
\description{Fitted values, residual extraction, simulation and print method for joint N-mixture models}
\usage{
\method{fitted}{jointNmix}(object, ...)
\method{print}{jointNmix}(x, round = TRUE, ...)
\method{residuals}{jointNmix}(object, type = c("ordinary", "standardized"), ...)
\method{simulate}{jointNmix}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, x}{fitted model object}
  \item{round}{logical. Round the estimates?}
  \item{type}{type of residuals to be returned. May be ordinary or standardized}
  \item{\dots}{not used}
}
\details{
\code{fitted.uniNmix} prints the fitted values for the model fit; \code{residuals.uniNmix} returns the ordinary (\deqn{r_i^o=y_i-\hat{\mu}_i}) or standardized (\deqn{r_i^s=\frac{y_i-\hat{\mu}_i}{\mbox{Var}(Y_i)}}) residuals; \code{simulate.uniNmix} simulates one sample from the fitted model.
}
\references{Moral, R.A., Hinde, J., \enc{Demétrio}{Demetrio}, C.G.B., Reigada, C. and Godoy, W.A.C. (submitted) Models for jointly estimating abundance of two unmarked site-associated species subject to imperfect detection.}
\author{
Rafael A. Moral <rafael_moral@yahoo.com.br>, Clarice G. B. \enc{Demétrio}{Demetrio} and John Hinde
}
\seealso{
\code{\link{Nmix}}
\code{\link{jointNmix}}
}