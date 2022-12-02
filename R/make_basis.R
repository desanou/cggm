# https://stackoverflow.com/questions/44438613/creating-canonical-basis-vectors-in-r

make_basis <- function(k, p = 10) {
  replace(numeric(p), k, 1)
  }
