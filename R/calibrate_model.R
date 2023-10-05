library(DEoptim)
library(adaptMCMC)
library(readxl)

source('R/likelihood.R')
source('R/prior_list.R')

# specify model type
# read in PBS job array index to specify type
aid <- Sys.getenv('PBS_ARRAYID')
# if PBS_ARRAYID doesn't exist, this should be passed as a command line argument or set interactively
if (aid == '') {
  if (!exists('type')) {
    args <- commandArgs(trailingOnly=TRUE)
    type <- args[1]
  }
} else {
  types <- c('random', 'mean', 'drift')
  id <- as.numeric(aid)
  type <- types[id]
}


# is this a pilot run or a production run?
pilot <- FALSE
if (pilot) {
  par_flag <- TRUE
  pilot_flag <- TRUE
} else {
  par_flag <- FALSE
  pilot_flag <- FALSE
}


if (par_flag) {library(parallel)}


# set list of parameters based on model structure
if (type == 'random') {
  parnames <- c('rho1', 'rho2', 'rho3', 'sigma_sq')
} else if (type == 'mean') {
  parnames <- c('eta', 'rho1', 'rho2', 'rho3', 'sigma_sq')
} else if (type == 'drift') {
  parnames <- c('a', 'b', 'rho1', 'rho2', 'rho3', 'sigma_sq')
}

# define prior distribution types and quantile
prior_df <- data.frame(name=parnames, type=character(length(parnames)), lower=numeric(length(parnames)), upper=numeric(length(parnames)), stringsAsFactors=FALSE)
for (p in parnames) {
  # get index
  i <- match(p, prior_df$name)
  # if the parameter is an AR coefficient
  if (grepl('rho', p)) {
    prior_df[i, 'type'] <- 'normal'
    prior_df[i, 'lower'] <- -3
    prior_df[i, 'upper'] <- 3
  } else if (p == 'sigma_sq') { # if the parameter is the variance
    prior_df[i, 'type'] <- 'log-normal'
    prior_df[i, 'lower'] <- 0.001
    prior_df[i, 'upper'] <- 1
  } else if (p == 'eta') { # if the parameter is the mean
    prior_df[i, 'type'] <- 'normal'
    prior_df[i, 'lower'] <- 2
    prior_df[i, 'upper'] <- 5
  } else if (p == 'a') { # if the parameter is the linear intercept
    prior_df[i, 'type'] <- 'normal'
    prior_df[i, 'lower'] <- 0
    prior_df[i, 'upper'] <- 3
  } else if (p == 'b') {# if the parameter is the linear slope
    prior_df[i, 'type'] <- 'normal'
    prior_df[i, 'lower'] <- -0.5
    prior_df[i, 'upper'] <- 0.5
  }
}
# convert prior data frame to list for use with log_pri()
priors <- create_prior_list(prior_df)

# load data
discount_xlsx <- read_excel("data/discount.xlsx")
dat <- c(read_excel("data/discount.xlsx", range="Real rates!K7:K231"))
# remove NAs
dat <- dat[!is.na(dat)]

# find MAP estimate as starting point for MCMC
if (!file.exists(paste0('output/map-', type, '.rds'))) {
  ## set DEoptim parameters for MAP estimate
  # set upper and lower bounds for each variable
  # we include all variables for ease of switching between models.
  all_parnames <- c('rho1', 'rho2', 'rho3', 'sigma_sq', 'eta', 'a', 'b')
  lbound <- c(-2, -2, -2, 0, 0, 0, -2)
  ubound <- c(2, 2, 2, 1, 7, 2, 0)
  
  # if DEoptim is run in parallel, set up cluster and evaluate
  # otherwise just evaluate
  # set DEoptim flags
  NP_scale <- 20 # number of particles per parameter
  trace <- FALSE
  n_iter <- 5000
  
  exportvars <- c('log_post', 'log_pri', 'log_lik_kf', 'compute_residuals', 'log_lik')
  if (par_flag) {
    map_out <- DEoptim::DEoptim(neg_log_post, lbound[match(parnames, all_parnames)], ubound[match(parnames, all_parnames)], control=list(parallelType=1, trace=trace, itermax=n_iter, NP=NP_scale*length(parnames), parVar=exportvars), parnames=parnames, priors=priors, dat=dat, type=type)
  } else {
    map_out <- DEoptim::DEoptim(neg_log_post, lbound[match(parnames, all_parnames)], ubound[match(parnames, all_parnames)], control=list(NP=NP_scale*length(parnames), itermax=n_iter, trace=trace), parnames=parnames, priors=priors, dat=dat, type=type)
  }
  saveRDS(map_out, paste0('output/map-', type, '.rds'))
} else {
  map_out <- readRDS(paste0('output/map-', type, '.rds'))
}

## set up MCMC run parameters
# set initial step size to be 10% of the standard deviation of the prior distributions when defined, or else some initial guess otherwise
stepsize <- numeric(length(parnames))
for (name in parnames) {
  if (priors[[name]][['type']] == 'uniform') {
    stepsize[match(name, parnames)] <- 0.1*(priors[[name]][['max']] - priors[[name]][['min']])/sqrt(12)
  } else if (priors[[name]][['type']] == 'normal') {
    stepsize[match(name, parnames)] <- 0.1*priors[[name]][['sd']]
  } else if (priors[[name]][['type']] == 'truncnorm') {
    stepsize[match(name, parnames)] <- 0.1*priors[[name]][['sd']]
  } else if (priors[[name]][['type']] == 'log-normal') {
    stepsize[match(name, parnames)] <- 0.1*sqrt((exp(priors[[name]][['sdlog']]^2)-1)*exp(2*priors[[name]][['meanlog']]+priors[[name]][['sdlog']]^2))
  }
}

init <- map_out$optim$bestmem
names(init) <- parnames

rate_accept_many <- 0.234 # optimal acceptance rate for multiple parameters
rate_accept_one <- 0.44 # optimal acceptance rate for one parameter
rate_accept <- rate_accept_many + (rate_accept_one - rate_accept_many)/length(parnames)

if (pilot_flag) {
  n_iter <- 5e4
  n_core <- 4
  
  adapt_start <- max(500, round(0.01*n_iter)) # when to start stepsize adaptation
  adapt_stop <- TRUE
  
  mcmc_out <- adaptMCMC::MCMC.parallel(log_post, n_iter, init=init, n.cpu=n_core, adapt=adapt_stop, n.start=adapt_start, acc.rate=rate_accept, gamma=0.51, list=TRUE, packages=c('truncnorm'), scale=stepsize, parnames=parnames, priors=priors, dat=dat, type=type)
  saveRDS(mcmc_out, paste0('output/mcmc-pilot-', type, '.rds'))
} else {
  n_iter <- 1e5
  adapt_start <- max(500, round(0.01*n_iter)) # when to start stepsize adaptation
  adapt_stop <- max(1e4, round(0.1*n_iter))
  
  mcmc_out <- adaptMCMC::MCMC(log_post, n_iter, init=init, adapt=adapt_stop, acc.rate=rate_accept, gamma=0.51, list=TRUE, scale=stepsize, parnames=parnames, priors=priors, dat=dat, type=type)
  saveRDS(mcmc_out, paste0('output/mcmc-', type, '.rds'))
}