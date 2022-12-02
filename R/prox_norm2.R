prox_norm2 <- function(lambda2_rho, vk) {

  nrm = norm(vk, type = "2")
  p = length(vk)

  if(nrm > lambda2_rho){
    z <- (1-lambda2_rho/nrm) * vk
  }else{
    z <- rep(0, p)
  }

  return(z)
}
