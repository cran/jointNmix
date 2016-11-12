residuals.jointNmix <-
function(object, type=c("ordinary","standardized"), ...) {
  type <- match.arg(type)
  f <- fitted.jointNmix(object)
  r1 <- as.numeric(object$sp1 - f$sp1)
  r2 <- as.numeric(object$sp2 - f$sp2)
  w1 <- as.numeric(1/sqrt(object$weights1))
  w2 <- as.numeric(1/sqrt(object$weights2))
  res <- switch(type, ordinary = list("sp1"=r1, "sp2"=r2),
                      standardized = list("sp1"=r1 * w1, "sp2"=r2 * w2))
  return(res)
}
