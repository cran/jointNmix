print.uniNmix <-
function(x, round=TRUE, ...) {
  mix <- x$mixture
  np1 <- ncol(x$Xp)
  nl1 <- ncol(x$Xl)
  nam <- names(coef(x))
  p.obj <- data.frame("Estimate"=x[[1]], "SE"=x[[2]], "z"=x[[3]], "p-value"=x[[4]])
  rownames(p.obj) <- nam
  if(round) p.obj <- round(p.obj, 4)
  p.obj[,4][p.obj[,4]<0.0001] <- "<0.0001"
  cat("Detection estimates", "\n")
  print(p.obj[1:np1,])
  cat("", "\n")
  cat("Abundance estimates", "\n")
  if(mix=="P") print(p.obj[(np1+1):(np1+nl1),])
  if(mix=="NB"|mix=="NeymanA") print(p.obj[(np1+1):(np1+nl1+1),])
}
