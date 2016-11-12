logLik.Nmix <-
function(object, ...) {
  llik <- -object$logLik
  ndf <- length(coef(object))
  nam.obj <- deparse(substitute(object))
  cat("logLik(", nam.obj, ") = ", round(llik, 4), " (df = ", ndf, ")", "\n", sep="")
  return(invisible(llik))
}
