anova.Nmix <-
function(object, object2, ...) {
  np1 <- length(coef(object))
  np2 <- length(coef(object2))
  llik1 <- logLik(object)
  llik2 <- logLik(object2)
  chi <- 2*abs(llik1 - llik2)
  df <- abs(np1 - np2)
  p <- pchisq(chi, df, lower.tail=FALSE)
  p <- ifelse(p<0.0001, "<0.0001", round(p, 4))
  ret <- data.frame("Chisq"=chi,"d.f."=df,"p-value"=p)
  return(ret)
}
