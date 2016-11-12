\name{plot-methods-Nmix}
\alias{plot.uniNmix}
\alias{plot.jointNmix}
\encoding{UTF-8}
\title{Plot Methods for Nmix objects}
\description{Plots residuals vs. fitted values or the conditional distributions of the abundances for each site}
\usage{
\method{plot}{uniNmix}(x, posterior = FALSE, layout, sites, restype, ...)
\method{plot}{jointNmix}(x, posterior = FALSE, layout, sites, restype, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{fitted model object}
  \item{posterior}{logical. If TRUE, plots the conditional distributions of the abundances. If FALSE, plots residuals vs. fitted values}
  \item{layout}{two-number vector indicating the partition of the graphical window in rows and columns, respectively}
  \item{sites}{index of sites to plot the conditional distribution of the abundances. If left unspecified, all sites are plotted}
  \item{restype}{type of residuals to be plotted. Defaults to "ordinary"}
  \item{\dots}{not used}
}
\references{Moral, R.A., Hinde, J., \enc{Demétrio}{Demetrio}, C.G.B., Reigada, C. and Godoy, W.A.C. (submitted) Models for jointly estimating abundance of two unmarked site-associated species subject to imperfect detection.}
\author{
Rafael A. Moral <rafael_moral@yahoo.com.br>, Clarice G. B. \enc{Demétrio}{Demetrio} and John Hinde
}
\seealso{
\code{\link{Nmix}}
\code{\link{jointNmix}}
}