source('R/filter.R')

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
  # run Kalman filter
  kf_out <- kf(r, S, Z, R, Q)
  # get length of data series
  N <- length(r)
  ll <- numeric(N)
  for (t in 1:N) {
    ll[t] <- log(abs(kf_out$G[,,t])) + t(kf_out$v[,,t]) %*% solve(kf_out$G[,,t]) %*% kf_out$v[,,t]
  }
  ld <- -0.5*(N*log(2*pi) + sum(ll)) # compute log-likelihood (up to a constant)
  # return list of log-likelihood and parameters needed for prediction
  list(log.density=ld, 
       prev.mean=kf_out$a[nrow(kf_out$a),], 
       prev.var=kf_out$P[,,dim(kf_out$P)[3]],
       S=S, Z=Z, R=R, Q=Q)
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
  # check the unit root constraint
  rho <- pars[grepl('rho*', parnames)]
  if (sum(rho) > 1) {
    return(-Inf)
  }
  # compute the log-prior density of the parameters
  lp <- log_pri(pars, parnames, priors)
  # if any of the parameters have zero prior density, return -Inf
  if (lp == -Inf) {
    return(-Inf)
  } 
  # otherwise compute the log-likelihood
  ll_out <- log_lik(pars, parnames, dat, type)
  # return the log-posterior density
  ll_out[['log.density']] + lp
}

neg_log_post <- function(pars, parnames, priors, dat, type) {
  -1*log_post(pars, parnames, priors, dat, type)
}