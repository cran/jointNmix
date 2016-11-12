print.jointNmix <-
function(x, round=TRUE, ...) {
  mix <- x$mixture
  includepsi <- x$includepsi
  np1 <- ncol(x$Xp1)
  np2 <- ncol(x$Xp2)
  nl1 <- ncol(x$Xl1)
  nl2 <- ncol(x$Xl2)
  npsi <- ncol(x$Xpsi)
  nam <- names(coef(x))
  p.obj <- data.frame("Estimate"=x[[1]], "SE"=x[[2]], "z"=x[[3]], "p-value"=x[[4]])
  rownames(p.obj) <- nam
  if(round) p.obj <- round(p.obj, 4)
  p.obj[,4][p.obj[,4]<0.0001] <- "<0.0001"
  cat("Detection estimates", "\n")
  cat("Species 1", "\n")
  print(p.obj[1:np1,])
  cat("Species 2", "\n")
  print(p.obj[(np1+1):(np1+np2),])
  cat("", "\n")
  cat("Abundance estimates", "\n")
  cat("Species 1", "\n")
  if(mix[1]=="P") print(p.obj[(np1+np2+1):(np1+np2+nl1),])
  if(mix[1]=="NB") print(p.obj[c((np1+np2+1):(np1+np2+nl1),(np1+np2+nl1+nl2+1)),])
  cat("Species 2", "\n")
  if(mix[2]=="P") print(p.obj[(np1+np2+nl1+1):(np1+np2+nl1+nl2),])
  if(mix[2]=="NB"&mix[1]=="P") print(p.obj[c((np1+np2+nl1+1):(np1+np2+nl1+nl2),(np1+np2+nl1+nl2+1)),])
  if(mix[2]=="NB"&mix[1]=="NB") print(p.obj[c((np1+np2+nl1+1):(np1+np2+nl1+nl2),(np1+np2+nl1+nl2+2)),])
  if(includepsi)  {
    cat("Species 2 psi estimates", "\n")
    print(p.obj[(nrow(p.obj)-npsi+1):(nrow(p.obj)),])
  }
}
