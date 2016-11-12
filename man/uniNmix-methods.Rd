\name{uniNmix-methods}
\alias{fitted.uniNmix}
\alias{print.uniNmix}
\alias{residuals.uniNmix}
\alias{simulate.uniNmix}
\encoding{UTF-8}
\title{Methods for univariate Nmix objects}
\description{Fitted values, residual extraction, simulation and print method for univariate N-mixture models}
\usage{
\method{fitted}{uniNmix}(object, ...)
\method{print}{uniNmix}(x, round = TRUE, ...)
\method{residuals}{uniNmix}(object, type = c("ordinary", "standardized"), ...)
\method{simulate}{uniNmix}(object, ...)
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
\author{
Rafael A. Moral <rafael_moral@yahoo.com.br>, Clarice G. B. \enc{Demétrio}{Demetrio} and John Hinde
}
\seealso{
\code{\link{Nmix}}
\code{\link{jointNmix}}
}