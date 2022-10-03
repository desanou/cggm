#' Initialize primal variables to identity matrices
init_primal <- function(size) {
  primal_var <- matrix(data = 0, nrow = size, ncol = size)
  diag(primal_var) <- 1
  return(primal_var)
}
