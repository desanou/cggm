admm_cggm <- function(Sigma, lambda = 0, rho1, rho2, sj, niter_max = 20, tol = 1e-2, w = 1){

  conv_admm <- FALSE

  Psi <- delta <- list()
  Q <- R <- list()
  U <- Z <- list()

  p <- ncol(Sigma)
  ndiff <- p*(p-1)/2 ## if no weigths

  I2 = matrix(0, 2, 2)
  diag(I2) = 1

  ## algo infos
  infos_primalTheta = infos_primalDelta = infos_primalPsi = NULL
  infos_dualZ = infos_dualU = NULL
  infos_theta = NULL
  infos_time = infos_conv_grad_descent = NULL

  ## Compute Q R
  ll=0
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      ll = ll+1
      mat <- matrix(0, p, p)
      diag(mat) = 1
      Q[[ll]] <- mat[-c(i,j), ]

      R[[ll]] <- t(mat[c(i,j), ])
    }
  }

  ## Init primal variables
  Theta = matrix(0, p, p)
  diag(Theta) = 1
  Theta_old = Theta

  for(l in 1:ndiff){
    Psi[[l]] <- Q[[l]] %*% Theta %*% R[[l]]
    dummy_Psi <-  Psi[[l]]
    delta[[l]] <- dummy_Psi[,1] - dummy_Psi[,2]
  }
  Psi_old = Psi
  delta_old = delta

  ## Init dual variables
  for(l in 1:ndiff){
    U[[l]] <- matrix(0, p-2, 2)
  }
  U_old = U

  for(l in 1:ndiff){
    Z[[l]] <- matrix(0, p-2, 1)
  }
  Z_old = Z


  iter_admm = 0

  while (!conv_admm && iter_admm <= niter_max) {

    time_iter <- system.time({
      iter_admm = iter_admm + 1
      conv_grad_descent <- FALSE
      iter_gd = 0

      ## (i) Update Θ
      while (!conv_grad_descent && iter_gd <= niter_max) {
        iter_gd = iter_gd + 1
        # Θ(k,j) ←Θ(k,j−1) − sj∇ΘL(Θ(k,j−1));

        sum_qtrp <- 0
        for (l in 1:length(Q)) {
          sum_qtrp = sum_qtrp + 2* t(Q[[l]]) %*% Q[[l]] %*% Theta %*% R[[l]] %*% t(R[[l]]) -
            2 * t(Q[[l]]) %*% Psi[[l]] %*% t(R[[l]]) +
            2 * t(Q[[l]]) %*% U[[l]] %*% t(R[[l]])
        }

        grad_lagrangian <- - solve(Theta) +
          Sigma +
          rho1/2 * sum_qtrp

        Theta_new <- Theta - sj * grad_lagrangian

        diff_theta <- norm(Theta_new - Theta, type = "F")

        cat("\n iter_grad: ", iter_gd, " diff_Theta: ", diff_theta)

        if(diff_theta <= tol){
          conv_grad_descent = TRUE
        }

        Theta = Theta_new
      }

      ## (ii) Update Ψl (∀l ∈M in parallel):
      rprimal_Psi = 0
      for (l in 1:length(Psi)) {
        Psi[[l]] <- (rho1 / (rho1 + 2 * rho2)) *
          (Q[[l]] %*% Theta %*% R[[l]] +
             U_old[[l]] +
             (rho2 / rho1) *
             (delta_old[[l]] - Z_old[[l]]) %*%
             t(make_basis(1,2) - make_basis(2, 2))
          ) %*%
          (I2 +
             (rho2 / rho1) * matrix(1, 2, 2))

        rprimal_Psi = rprimal_Psi + norm(Psi[[l]] - Psi_old[[l]], type = 'F')
      }

      ## (iii) Update δl (∀l ∈M in parallel): δ(k)
      rprimal_delta = 0
      for (l in 1:length(delta)) {
        Dvec <- Psi[[l]][,1] - Psi[[l]][,2]; dim(Dvec) <- dim(Z_old[[l]])
        delta[[l]] <- prox_norm2(lambda * w / rho2, Dvec + Z_old[[l]])

        rprimal_delta = rprimal_delta + norm(delta[[l]] - delta_old[[l]], type = "2")
      }

      ## (iv) Update Ul (∀l ∈M in parallel):
      rdual_U = 0
      for (l in 1:length(U)) {
        U[[l]] <- U_old[[l]] + (Q[[l]] %*% Theta %*% R[[l]] - Psi[[l]])

        rdual_U = rdual_U + norm(U[[l]] - U_old[[l]], type = 'F')
      }

      ## (v) Update zl (∀l ∈M in parallel):
      rdual_Z = 0
      for (l in 1:length(Z)) {
        Dvec <- Psi[[l]][,1] - Psi[[l]][,2]; dim(Dvec) <- dim(Z_old[[l]])
        Z[[l]] <- Z_old[[l]] + Dvec - delta[[l]]

        # cat("\n iter_admm: ", iter_admm, " Z", l, ": ", Z[[l]])

        rdual_Z = rdual_Z + norm(Z[[l]] - Z_old[[l]], type = "2")
      }
    })

    rprimal_Theta <- norm(Theta - Theta_old, type = 'F')

    cat("\n iter_admm: ", iter_admm,
        " primal res Theta: ", rprimal_Theta,
        " Psi: ", rprimal_Psi,
        " delta: ", rprimal_delta,
        " dual res U: ", rdual_U,
        " Z: ", rdual_Z)

    Theta_old = Theta
    Psi_old = Psi
    delta_old = delta
    Z_old = Z
    U_old = U

    infos_theta[[iter_admm]] = Theta
    infos_cost = NULL
    infos_primalTheta[[iter_admm]] <- rprimal_Theta
    infos_primalDelta[[iter_admm]] <- rprimal_delta
    infos_primalPsi[[iter_admm]] <- rprimal_Psi
    infos_dualU[[iter_admm]] <- rdual_U
    infos_dualZ[[iter_admm]] <- rdual_Z
    infos_time[[iter_admm]] <- time_iter
    infos_conv_grad_descent[[iter_admm]] <- conv_grad_descent

    if (max(rprimal_Theta, rprimal_delta, rprimal_Psi, rdual_U, rdual_Z) <= tol){
      conv_admm = TRUE
      break
    }
  }

  res <- list(conv_admm = conv_admm,
              iter_admm = iter_admm,
              infos_theta = infos_theta,
              infos_cost = NULL,
              infos_dualU = infos_dualU,
              infos_dualZ = infos_dualZ,
              infos_primalPsi = infos_primalPsi,
              infos_primalDelta = infos_primalDelta,
              infos_time = infos_time,
              infos_primalTheta = infos_primalTheta,
              infos_conv_grad_descent = infos_conv_grad_descent,
              Theta = Theta,
              Psi = Psi,
              delta = delta)

  return(res)
}
