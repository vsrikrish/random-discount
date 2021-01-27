## set up prior list based on passed in data frame

create_prior_list <- function(prior_df) {
  parnames <- prior_df[, 'name']
  priors <- vector('list', length(parnames))
  names(priors) <- parnames
  
  for (i in 1:nrow(prior_df)) {
    name <- prior_df[i, 'name']
    priors[[name]] <- list(type=prior_df[i, 'type'])
    if (priors[[name]][['type']] == 'uniform') {
      priors[[name]][['dens.fun']] <- 'dunif'
      priors[[name]][['min']] <- prior_df[i, 'lower']
      priors[[name]][['max']] <- prior_df[i, 'upper']
    } else if (priors[[name]][['type']] == 'normal') {
      priors[[name]][['dens.fun']] <- 'dnorm'
      priors[[name]][['mean']] <- mean(c(prior_df[i, 'lower'], prior_df[i, 'upper']))
      priors[[name]][['sd']] <- (prior_df[i, 'upper'] - prior_df[i, 'lower'])/(qnorm(0.975) - qnorm(0.025))
    } else if (priors[[name]][['type']] == 'log-normal') {
      priors[[name]][['dens.fun']] <- 'dlnorm'
      priors[[name]][['meanlog']] <- mean(c(log(prior_df[i, 'lower']), log(prior_df[i, 'upper'])))
      priors[[name]][['sdlog']] <- (log(prior_df[i, 'upper']) - log(prior_df[i, 'lower'])) / (qnorm(0.975) - qnorm(0.025))
    } else if (priors[[name]][['type']] == 'truncnorm') {
      require(truncnorm)
      priors[[name]][['dens.fun']] <- 'dtruncnorm'
      priors[[name]][['a']] <- 0
      priors[[name]][['b']] <- Inf
      priors[[name]][['mean']] <- mean(c(prior_df[i, 'lower'], prior_df[i, 'upper']))
      priors[[name]][['sd']] <- (prior_df[i, 'upper'] - prior_df[i, 'lower'])/(qnorm(0.975) - qnorm(0.025))
    }
  }
  
  priors
}