plot.uniNmix <-
function(x, posterior=FALSE, layout, sites, restype, ...) {
  if(!posterior) {
    if(missing(restype)) restype <- "ordinary"
    fv <- as.numeric(fitted.uniNmix(x)[[1]])
    res <- residuals.uniNmix(x, type=restype)
    plot(fv, res, xlab="Fitted values", ylab="Residuals")
    abline(h=0, lty=2)
    lines(lowess(fv, res), col=2)
  } else {
    rf <- getranef.uniNmix(x, distr=TRUE)[[1]]
    colnames(rf) <- 1:ncol(rf)
    rf <- rf[,apply(rf, 2, sum)>1e-6]
    if(missing(sites)) sites <- 1:nrow(rf)
    if(missing(layout)) layout <- c(4, ceiling(nrow(rf)/4))
    par(mfrow=layout, mar=c(2,2,2,2))
    for(i in sites) plot(rf[i,], ty="h", ylim=c(0, max(rf)), main=eval(i))
    par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
  }
}
