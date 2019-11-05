## log_pri
# compute prior densities
log_pri <- function(pars, parnames, priors) {
  
  # this function evaluates the log-prior density for a given parameter
  log_dens <- function(p) {
    val <- pars[match(p, parnames)]
    do.call(match.fun(priors[[p]][['dens.fun']]),
            c(list(x=val, log=TRUE),
              priors[[p]][-which(names(priors[[p]]) %in% c('type', 'dens.fun'))])
    )
  }
  # evaluate log-prior densities for each parameter
  lp <- vapply(parnames, log_dens, numeric(1))
  
  # return sum of log-priors
  sum(lp)
}

## compute_residuals
# compute the residuals of the data given the parameters and the model type
compute_residuals <- function(pars, parnames, dat, type) {
  if (type == 'random') {
    dat - head(dat, 1)
  } else if (type == 'mean') {
    dat - log(pars[match('eta', parnames)])
  } else if (type == 'drift') {
    N <- length(dat)
    dat - (pars[match('a', parnames)] + pars[match('b', parnames)] * 1:N)
  }
}

# run the Kalman filter given the passed matrices and return the log-likelihood
log_lik_kf <- function(r, S, Z, R, Q) {
  # get length of data series
  N <- length(r)
  # find initial covariance and state estimate
  P0 <- matrix(solve(diag(dim(S)[1]^2) - kronecker(S, S)) %*% as.numeric(R %*% Q %*% t(R)), ncol=3)
  a0 <- c(0, 0, 0)
  
  # set up storage for estimates
  v <- array(NA, dim=c(1, 1, N))
  G <- array(NA, dim=c(1, 1, N))
  P <- array(NA, dim=c(3, 3, N))
  a <- matrix(0, nrow=N, ncol=3)
  ll <- numeric(N)
  
  # store initial values for the filter
  a[1,] <- a0
  P[,,1] <- P0

  # run the filter and compute log-likelihod components
  for (t in 2:N) {
    v[,,t-1] <- r[t-1] - Z %*% a[t-1,]
    G[,,t-1] <- Z %*% P[,,t-1] %*% t(Z)
    ll[t-1] <- log(abs(G[,,t-1])) + t(v[,,t-1]) %*% solve(G[,,t-1]) %*% v[,,t-1]
    at <- a[t-1,] + (P[,,t-1] %*% t(Z)  %*% solve(G[,,t-1]) %*% v[,,t-1])
    Pt <- P[,,t-1] - (P[,,t-1] %*% t(Z) %*% solve(G[,,t-1]) %*% Z %*% P[,,t-1])
#    ll[t-1] <- log(abs(G)) + (v %*% solve(G) %*% v) # compute the log 
    a[t,] <- S %*% at
    P[,,t] <- (S %*% Pt %*% solve(S)) + (R %*% Q %*% t(R))
  }
  # compute log-likelihood component from last time step
  v[,,N] <- r[t-1] - Z %*% a[t-1,]
  G[,,N] <- Z %*% P[,,t-1] %*% t(Z) 
  ll[N] <- log(abs(G[,,N])) + t(v[,,N]) %*% solve(G[,,N]) %*% v[,,N]
  -0.5*(N*log(2*pi) + sum(ll)) # return log-likelihood (up to a constant)
}

## compute the log-likelihood of the parameters for the given model type
log_lik <- function(pars, parnames, dat, type) {
  r <- compute_residuals(pars, parnames, dat, type) # compute residuals
  # set up state space matrices from sampled parameters
  rho <- pars[grepl('rho*', parnames)]
  S <- as.matrix(cbind(rho, rbind(diag(2), numeric(2)))) #companion matrix
  Z <- matrix(0, nrow=1, ncol=3) # observation matrix
  Z[1,1] <- 1
  R <- matrix(0, nrow=3, ncol=1)
  R[1,1] <- 1
  Q <- diag(pars[match('sigma_sq', parnames)], 1)
  # compute log-likelihood and return
  log_lik_kf(r, S, Z, R, Q)
} 

log_post <- function(pars, parnames, priors, dat, type) {
  # compute the log-prior density of the parameters
  lp <- log_pri(pars, parnames, priors)
  # if any of the parameters have zero prior density, return -Inf
  if (lp == -Inf) {
    return(-Inf)
  } 
  # otherwise compute the log-likelihood
  ll <- log_lik(pars, parnames, dat, type)
  # return the log-posterior density
  ll + lp
}

neg_log_post <- function(pars, parnames, priors, dat, type) {
  -1*log_post(pars, parnames, priors, dat, type)
}