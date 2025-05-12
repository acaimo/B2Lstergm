#' Title Simulate Network Statistics from a Separable 2-layer TERGM
#'
#' @param formula 2-layer ERGM model formula
#' @param y_t list of multi-layer networks
#' @param L.x_F list of observed x_F networks
#' @param L.x_P list of observed x_P networks
#' @param L.z_F initial z_F state (multi-layer network)
#' @param L.z_P initial z_P state (multi-layer network)
#' @param times number of time points, including t_0
#' @param model.dim model dimension, including constraint: L(~edges, ~(active))
#' @param coef_F parameter zeta_F without constraint: (Inf,...
#' @param coef_P parameter zeta_P without constraint: (Inf,...
#' @param ergm_sim_iters number of iterations for ERGM simulation
#' @param network logical, if TRUE saves the simulated networks z_0:T
#' @param network.stats logical, if TRUE saves the simulated network statistics of z_0:T
#' @param ... additional arguments
#'
#' @returns list of network statistics and network objects
#' @import ergm ergm.multi mvtnorm
#' @importFrom stats cov runif simulate
#' @export
simulate_stats <- function(formula,
                           y_t = y_t,
                           L.x_F = L.x_F,
                           L.x_P = L.x_P,
                           L.z_F = L.z_F,
                           L.z_P = L.z_P,
                           times = NULL,
                           model.dim = NULL,
                           coef_F = NULL,
                           coef_P = NULL,
                           ergm_sim_iters = 1000,
                           network = FALSE,
                           network.stats = FALSE, ...) {

  # set up output storage
  sz._F <- sz._P <- matrix(NA, model.dim, times)
  sz. <- if (network.stats) matrix(NA, model.dim, times) else NULL

  for(t in 2:times){
    #
    # Simulation
    #
    suppressMessages(
      y._F <- simulate(formula,
                       coef = c(Inf, coef_F), # Formation
                       nsim = 1,
                       constraints = ~Dyads(vary = ~dyadcov(Layer(L.x_F[[t - 1]], L.x_F[[t - 1]]))),
                       basis = Layer(positive = L.z_F, active = L.x_F[[t - 1]]),
                       control = control.simulate(MCMC.burnin = ergm_sim_iters),
                       output = 'network')
    )
    suppressMessages(
      y._P <- simulate(formula,
                       coef = c(Inf, coef_P), # Persistence
                       nsim = 1,
                       constraints = ~Dyads(vary = ~dyadcov(Layer(L.x_P[[t - 1]], L.x_P[[t - 1]]))),
                       basis = Layer(positive = L.z_P, active = L.x_P[[t - 1]]),
                       control = control.simulate(MCMC.burnin = ergm_sim_iters),
                       output = 'network')
    )
    sz._F[, t] <- summary(formula, basis = y._F)
    sz._P[, t] <- summary(formula, basis = y._P)
    #
    # Update y_t[[t]], L.z_F and L.z_P
    #
    if(t < times) {
      # Update only if there are exogenous dynamic network statistics
      # y_t[[t]][, ] <- y._P[, ] | (y._F[, ] - y_t[[t - 1]][, ])
      L.z_F <- uncombine_network(y._F)[[1]]
      L.z_P <- uncombine_network(y._P)[[1]]
    }
    if(network.stats) {
      y_t[[t]][, ] <- y._P[, ] | (y._F[, ] - y_t[[t - 1]][, ])
      sz.[, t] <- summary(formula, basis = y_t[[t]])
    }
    if(network) y_t[[t]][, ] <- y._P[, ] | (y._F[, ] - y_t[[t - 1]][, ])
  }
  list(# including L(~edges, ~(active))
       # excluding t_0
       sz._F = unname(sz._F[, -1]),
       sz._P = unname(sz._P[, -1]),
       y_t = if(network) y_t else NULL,
       sz. = if (network.stats) sz.[, -1] else NULL)
}


