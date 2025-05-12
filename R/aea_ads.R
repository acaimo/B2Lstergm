#' Title Approximate Exchange Algorithm (AEA) with Adaptive Direction Sampling (ADS)
#'
#' @param formula 2-layer ERGM model formula
#' @param y_t list of observed multi-layer networks
#' @param L.x_F list of observed x_F networks
#' @param L.x_P list of observed x_P networks
#' @param L.z_F list of observed z_F networks
#' @param L.z_P list of observed z_P networks
#' @param sy observed statistics, i.e., list(sz$F_obs, sz$P_obs), include the L(~edges, ~(active))
#' @param iters overall number of iterations for each MCMC chain
#' @param burn.in number of burnin iterations for each MCMC chain
#' @param ergm_sim_iters number of iterations for ERGM simulation
#' @param mu.prior_F prior mean for zeta_F
#' @param Sigma.prior_F prior covariance for zeta_F
#' @param mu.prior_P prior mean for zeta_P
#' @param Sigma.prior_P prior covariance for zeta_P
#' @param n_chains number of MCMC chains
#' @param gamma_F ADS move factor for the F process
#' @param gamma_P ADS move factor for the P process
#' @param sigma.prop_F diagonal entry for the multivariate Normal proposal for the F process
#' @param sigma.prop_P diagonal entry for the multivariate Normal proposal for the P process
#'
#' @returns posterior samples of theta_F and theta_P
#' @import ergm ergm.multi mvtnorm
#' @importFrom stats cov runif simulate
#' @export
aea_ads <- function(formula,
                    y_t = y_t,
                    L.x_F = L.x_F,
                    L.x_P = L.x_P,
                    L.z_F = L.z_F,
                    L.z_P = L.z_P,
                    sy = sy,

                    iters = 50000,
                    burn.in = 100,
                    ergm_sim_iters = 1000,
                    # Prior specification
                    mu.prior_F = NULL,
                    Sigma.prior_F = NULL,
                    mu.prior_P = NULL,
                    Sigma.prior_P = NULL,
                    # ADS setting
                    n_chains = 10,
                    gamma_F = 0.5,
                    gamma_P = 0.5,
                    sigma.prop_F = 0.025,
                    sigma.prop_P = 0.025){
  #
  model.dim <- dim(sy$F_obs)[1] - 1 # remove first term L(~edges, ~(active))
  times <- length(L.x_F) + 1        # including t_0
  #
  if(is.null(mu.prior_F))    mu.prior_F <- rep(0, model.dim)
  if(is.null(Sigma.prior_F)) Sigma.prior_F <- diag(10, model.dim)
  if(is.null(mu.prior_P))    mu.prior_P <- rep(0, model.dim)
  if(is.null(Sigma.prior_P)) Sigma.prior_P <- diag(10, model.dim)
  #
  Sigma.prop_F <- diag(sigma.prop_F, model.dim)
  Sigma.prop_P <- diag(sigma.prop_P, model.dim)
  #
  theta_F <- array(NA, c(iters - burn.in, model.dim, n_chains))
  theta_P <- array(NA, c(iters - burn.in, model.dim, n_chains))
  #
  Delta_F <- matrix(0, length(seq(1001, 50000, 100)), model.dim)
  Delta_P <- matrix(0, length(seq(1001, 50000, 100)), model.dim)
  ii <-0
  #
  coef_F <- matrix(mu.prior_F +
                     runif(model.dim * n_chains,
                           min = -0.1, max = 0.1),
                   model.dim, n_chains)
  #
  coef_P <- matrix(mu.prior_P +
                     runif(model.dim * n_chains,
                           min = -0.1, max = 0.1),
                   model.dim, n_chains)
  #
  message(" > MCMC start")
  for (i in 1:iters) {
    for (h in 1:n_chains) {
      #
      coef._F <- coef_F[, h] +
        gamma_F * diff(c(coef_F)[sample(seq(1, n_chains)[-h], 2)]) +
        rmvnorm(1, sigma = Sigma.prop_F)[1, ]
      #
      coef._P <- coef_P[, h] +
        gamma_P * diff(c(coef_P)[sample(seq(1, n_chains)[-h], 2)]) +
        rmvnorm(1, sigma = Sigma.prop_P)[1, ]
      #
      sy. <- simulate_stats(formula,
                            y_t = y_t,
                            L.x_F = L.x_F,
                            L.x_P = L.x_P,
                            L.z_F = L.z_F,
                            L.z_P = L.z_P,
                            times = times,
                            model.dim = model.dim + 1,
                            coef_F = coef._F,
                            coef_P = coef._P,
                            ergm_sim_iters = ergm_sim_iters)
      #
      # NB: [-1] => remove offset parameter
      delta_F <- rowSums(sy.$sz._F - sy$F_obs)[-1]
      delta_P <- rowSums(sy.$sz._P - sy$P_obs)[-1]
      #
      prior_F <- dmvnorm(rbind(coef._F, coef_F[, h]),
                         mean = mu.prior_F,
                         sigma = Sigma.prior_F,
                         log = TRUE)
      #
      prior_P <- dmvnorm(rbind(coef._P, coef_P[, h]),
                         mean = mu.prior_P,
                         sigma = Sigma.prior_P,
                         log = TRUE)
      #
      beta_F <- (coef_F[, h] - coef._F) %*% delta_F + prior_F[1] - prior_F[2]
      beta_P <- (coef_P[, h] - coef._P) %*% delta_P + prior_P[1] - prior_P[2]
      #
      if (c(beta_F) >= log(runif(1))) {
        coef_F[, h] <- coef._F
      }
      if (c(beta_P) >= log(runif(1))) {
        coef_P[, h] <- coef._P
      }
    }
    #
    if (i > burn.in){
      theta_F[i - burn.in, , ] <- coef_F
      theta_P[i - burn.in, , ] <- coef_P
    }
    #
    if (i %in% seq(100, 50000, 100)){
      message(" >", i)
    }
  }
  list(theta_F = theta_F, theta_P = theta_P)
}
