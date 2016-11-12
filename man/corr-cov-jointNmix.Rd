\name{corr-cov-jointNmix}
\alias{getcorr.jointNmix}
\alias{getcov.jointNmix}
\encoding{UTF-8}
\title{Extract correlations and covariances}
\description{Extract estimated correlations and covariances between observed abundances for joint N-mixture models}
\usage{
getcorr.jointNmix(obj)
getcov.jointNmix(obj)
}
\arguments{
  \item{obj}{fitted model object}
}
\details{
The correlations under each specification of the latent abundance distributions are displayed in Moral et al. (submitted)'s Table 1. The covariances are derived as Supplementary information.
}
\references{Moral, R.A., Hinde, J., \enc{Demétrio}{Demetrio}, C.G.B., Reigada, C. and Godoy, W.A.C. (submitted) Models for jointly estimating abundance of two unmarked site-associated species subject to imperfect detection.}
\author{
Rafael A. Moral <rafael_moral@yahoo.com.br>, Clarice G. B. \enc{Demétrio}{Demetrio} and John Hinde
}

\seealso{
\code{\link{jointNmix}}
}