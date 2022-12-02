library(Matrix)
library(mvtnorm)

n = 10
K = 3
p = 9
rho = 0.85
blocs <- list()
for (j in 1:K) {
  bloc <- matrix(rho, nrow = p/K, ncol = p/K)
  for(i in 1:(p/K)) { bloc[i,i] <- 1 }
  blocs[[j]] <- bloc
}

mat.covariance <- bdiag(blocs)
mat.covariance
## Sparsity
## Sim
set.seed(2021)
X <- rmvnorm(n, mean = rep(0,p), sigma = as.matrix(mat.covariance))
X <- scale(X)

Sigma_hat = cov(X)

Omega_cggm <- admm_cggm(Sigma = Sigma_hat,
                        rho1 = 0.1,
                        rho2 = 0.1,
                        sj = 0.1,
                        lambda = 0.01,
                        niter_max = 1e5,
                        tol = 1e-2)

find_clusters(Omega_cggm$delta, p, thresh_fuse = 1)
