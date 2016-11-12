\name{Nmix-methods}
\alias{AIC.Nmix}
\alias{anova.Nmix}
\alias{coef.Nmix}
\alias{logLik.Nmix}
\encoding{UTF-8}
\title{Methods for Nmix objects}
\description{AIC, anova, coef, and logLik methods for univariate and joint N-mixture models}
\usage{
\method{AIC}{Nmix}(object, ...)
\method{anova}{Nmix}(object, object2, ...)
\method{coef}{Nmix}(object, ...)
\method{logLik}{Nmix}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{fitted model object}
  \item{object2}{second model object for the \code{anova} method}
  \item{\dots}{not used}
}
\details{
\code{AIC.Nmix} provides the Akaike Information Criterion for the model fit; \code{anova.Nmix} performs a likelihood-ratio test between two nested model fits; \code{coef.Nmix} prints the estimates of the model fit; \code{logLik.Nmix} prints the log-likelihood of the model fit.
}
\author{
Rafael A. Moral <rafael_moral@yahoo.com.br>, Clarice G. B. \enc{Demétrio}{Demetrio} and John Hinde
}
\seealso{
\code{\link{Nmix}}
\code{\link{jointNmix}}
}