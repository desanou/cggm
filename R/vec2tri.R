vec2tri <-function(k,p) {
  i <- ceiling(0.5*(2*p-1 - sqrt((2*p-1)^2 - 8*k)))
  j <- k - p*(i-1) + i*(i-1)/2 + i
  return(as.matrix(cbind(i,j)))
}
