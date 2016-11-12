dneymanA <- Vectorize(function(x, lambda1, lambda2, K, log=FALSE) {
  j <- 0:K
  dna <- exp(-lambda1)*lambda2^x/factorial(x)*sum(lambda1*exp(-lambda2)^j*j^x/(factorial(j)))
  if(log) return(log(dna)) else return(dna)
}, "x")
