\name{dneymanA}
\alias{dneymanA}
\encoding{UTF-8}
\title{The Neyman-A probability function}
\description{Computes the probability function of the Neyman-A distribution}
\usage{
dneymanA(x, lambda1, lambda2, K, log = FALSE)
}
\arguments{
  \item{x}{vector of values}
  \item{lambda1, lambda2}{parameters of the distribution}
  \item{K}{truncation value for the infinite summation}
  \item{log}{logical. If TRUE, the logarithm of the probabilities is returned}
}
\details{
The Neyman-A distribution has probability function \deqn{\frac{e^{-\lambda_1}\lambda_2^{x}}{x!}\sum_{k=0}^\infty\frac{(\lambda_1e^{-\lambda_2})^kk^x}{k!}} and is an overdispersion model. The summation is truncated to K.
}
\author{
Rafael A. Moral <rafael_moral@yahoo.com.br>, Clarice G. B. \enc{Demétrio}{Demetrio} and John Hinde
}\examples{
x <- 0:10
dneymanA(x, lambda1 = 2, lambda2 = 1.5, K = 50)
}