#####################################################
# random_discount.R                                 #
# This file includes functions to implement random  #
#   discounting, as in Newell and Pizer (2007).     #
# We only use the best-estimate of parameters, as   #
#   we did not have access to the data or the joint #
#   covariance matrix of the parameters.            #
# Functions are provided for both mean-reverting    #
#   and random-walk models without time trends.     #
# We initialize the models by assuming a 4% rate    #
#   of return or zero innovations by default in     #
#   2017-2018.                                      #
#####################################################

############## drift_discount #####################
## discount rate with linear trend
##    Inputs:
##            1) d (data time series, on a log scale)
##            2) N (integer length of desired series, starting in 2019)
##            3) n (integer number of replicates)
##            4) pars (named vector of parameter values)
####################################################

drift_discount <- function(d, N=100, pars=c(a=1.9349, b=-0.0058, rho1=1.6965, rho2=-0.9754, rho3=0.2388, sigma_sq=0.0034), n=nrow(pars), parnames=names(pars)) {
  
 pars <- as.matrix(pars) 
 colnames(pars) <- parnames
 offset <- length(d)-1
 # compute the trend
 tr <- sweep(outer(pars[,'b'], (offset-2):(N + offset), FUN='*'), 1, pars[,'a'], '+')

 # initialize storage for innovations
 eps <- matrix(NA, ncol=N+3, nrow=n)
 # add initial values
 eps[,1:3] <- matrix(tail(d, 3) - tr[1:3], nrow=n, ncol=3, byrow=TRUE)
 # for each subsequent year, use model to determine innovations
 for (i in 4:ncol(eps)) {
   eps[,i] <- pars[,'rho1'] * eps[,i-1] + pars[,'rho2'] * eps[,i-2] + pars[,'rho3'] * eps[,i-3] + rnorm(n=n, sd=sqrt(pars[,'sigma_sq']))
 }
 # add the trend to the innovation process and exponentiate
 exp(eps[,4:ncol(eps)] + tr[,4:ncol(tr)])
}


############## mean_discount #####################
## mean-reverting discount rate
##    Inputs:
##            1) d (data time series, on a log scale)
##            2) N (integer length of desired series, starting in 2019)
##            3) n (integer number of replicates)
##            4) pars (named vector of parameter values)
####################################################
mean_discount <- function(d,
                         N = 100,
                         pars = c(eta = 3.4,
                                  rho1 = 1.7371,
                                  rho2 = -1.0270,
                                  rho3 = 0.2806,
                                  sigma_sq = 0.0033),
                         n = nrow(pars),
                         parnames=names(pars)
){
  
  pars <- as.matrix(pars) 
  colnames(pars) <- parnames
  
  # initialize storage for innovations
  eps <- matrix(NA, ncol=N+3, nrow=n)
  
  # add initial values
  eps[,1:3] <- matrix(t(apply(pars, 1, function(p) tail(d, 3) - log(p['eta']))), nrow=n, ncol=3, byrow=TRUE)
  # for each subsequent year, use mean-reverting model to determine discount rates
  for (i in 4:ncol(eps)) {
    eps[,i] <- pars[,'rho1'] * eps[,i-1] +
      pars[,'rho2'] * eps[,i-2] +
      pars[,'rho3'] * eps[,i-3] +
      rnorm(n=n, sd=sqrt(pars[,'sigma_sq']))
  }
  
  # compute the bond rates starting at the mean value
  exp(sweep(eps[, 4:ncol(eps)], 1, log(pars[,'eta']), FUN='+'))
}

############## random_discount #####################
## random-walk discount rate
##    Inputs:
##            1) N (integer length of desired bond rate series, starting in 2019)
##            2) n (integer number of replicates)
##            3) pars (named vector of parameter values)
####################################################
random_discount <- function(d,
                        N = 100,
                        n = 1e5,
                        pars = c(rho1 = 1.7429,
                                 rho2 = -1.0455,
                                 rho3 = 0.3010,
                                 sigma_sq = 0.0033),
                        parnames=names(pars)
){
  pars <- as.matrix(pars) 
  colnames(pars) <- parnames
  
  # initialize storage for innovations
  eps <- matrix(NA, ncol=N+3, nrow=n)
  # add initial values based on deviation from initial value
  eps[,1:3] <- matrix(tail(d, 3) - head(d[!is.na(d)], 1), nrow=n, ncol=3, byrow=TRUE)
  # for each subsequent year, use mean-reverting model to determine discount rates
  for (i in 4:ncol(eps)) {
    eps[,i] <- pars[,'rho1'] * eps[,i-1] +
      pars[,'rho2'] * eps[,i-2] +
      pars[,'rho3'] * eps[,i-3] +
      rnorm(n=n, sd=sqrt(pars[,'sigma_sq']))
  }
  
  # compute the bond rates starting at the mean value
  exp(head(d[!is.na(d)],1) + eps[,4:ncol(eps)])
}

############## compute_factor #####################
## compute discount factors
##    Inputs:
##            1) dr (matrix of discount rates. Rows should correspond to a given realization).
####################################################
compute_factor <- function(dr, t=NULL) {
 # compute the discount factor
 if (is.matrix(dr)) {
   df <- cbind(1, t(exp(-1*apply(dr / 100, 1, cumsum))))
 } else {
   df <- exp(-dr*seq(0, t, 1))
 }
 df
}

############## ce_discount.R #####################
## certain-equivalent discount factor
##    Inputs:
##            1) df (matrix of discount factors. Rows should correspond to a given realization).
####################################################
ce_discount <- function(df) {
  # compute the certainty-equivalent discount factor, which is the expectation
  colMeans(df)
}
