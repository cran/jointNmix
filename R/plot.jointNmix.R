plot.jointNmix <-
function(x, posterior=FALSE, layout, sites, restype, ...) {
  if(!posterior) {
    if(missing(restype)) restype <- "ordinary"
    fv1 <- as.numeric(fitted.jointNmix(x)[[1]])
    res1 <- as.numeric(residuals.jointNmix(x, type=restype)[[1]])
    fv2 <- as.numeric(fitted.jointNmix(x)[[2]])
    res2 <- as.numeric(residuals.jointNmix(x, type=restype)[[2]])
    par(mfrow=c(1,2))
    plot(fv1, res1, xlab="Fitted values", ylab="Residuals", main="Species 1")
    abline(h=0, lty=2)
    lines(lowess(fv1, res1), col=2)
    plot(fv2, res2, xlab="Fitted values", ylab="Residuals", main="Species 2")
    abline(h=0, lty=2)
    lines(lowess(fv2, res2), col=2)
  } else {
    rf <- getranef.jointNmix(x, distr=TRUE)
    rf1 <- rf[[1]]
    colnames(rf1) <- 1:ncol(rf1)
    rf1 <- rf1[,apply(rf1, 2, sum)>1e-6]
    rf2 <- rf[[2]]
    colnames(rf2) <- 1:ncol(rf2)
    rf2 <- rf2[,apply(rf2, 2, sum)>1e-6]
    if(missing(sites)) sites <- 1:nrow(rf1)
    if(missing(layout)) layout <- c(4, ceiling(nrow(rf1)/4))
    par(mfrow=layout, mar=c(2,2,2,2))
    for(i in sites) plot(rf1[i,], type="h", ylim=c(0, max(rf1)), main=eval(i))
    par(mfrow=layout, mar=c(2,2,2,2), ask=TRUE)
    for(i in sites) plot(rf2[i,], type="h", ylim=c(0, max(rf2)), main=eval(i))
    par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), ask=F)
  }
}
