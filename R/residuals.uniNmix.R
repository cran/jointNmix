residuals.uniNmix <-
function(object, type=c("ordinary","standardized"), ...) {
  type <- match.arg(type)
  f <- fitted.jointNmix(object)
  r1 <- as.numeric(object$sp1 - f$sp1)
  w1 <- as.numeric(1/sqrt(object$weights1))
  res <- switch(type, ordinary = r1,
                standardized = r1 * w1)
  return(res)
}
