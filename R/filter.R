# run the Kalman filter given the passed matrices and data
kf <- function(r, S, Z, R, Q) {
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
  Pt <- array(NA, dim=c(3, 3, N))
  at <- matrix(0, nrow=N, ncol=3)
  
  
  # store initial values for the filter
  a[1,] <- a0
  P[,,1] <- P0
  
  # run the filter and compute log-likelihod components
  for (t in 2:N) {
    v[,,t-1] <- r[t-1] - Z %*% a[t-1,]
    G[,,t-1] <- Z %*% P[,,t-1] %*% t(Z)
    at[t-1,] <- a[t-1,] + (P[,,t-1] %*% t(Z)  %*% solve(G[,,t-1]) %*% v[,,t-1])
    Pt[,,t-1] <- P[,,t-1] - (P[,,t-1] %*% t(Z) %*% solve(G[,,t-1]) %*% Z %*% P[,,t-1])
    #    ll[t-1] <- log(abs(G)) + (v %*% solve(G) %*% v) # compute the log 
    a[t,] <- S %*% at[t-1,]
    P[,,t] <- (S %*% Pt[,,t-1] %*% solve(S)) + (R %*% Q %*% t(R))
  }
  # compute log-likelihood component from last time step
  v[,,N] <- r[N] - Z %*% a[N,]
  G[,,N] <- Z %*% P[,,N] %*% t(Z) 
  at[N,] <- a[N,] + (P[,,N] %*% t(Z)  %*% solve(G[,,N]) %*% v[,,N])
  Pt[,,N] <- P[,,N] - (P[,,t] %*% t(Z) %*% solve(G[,,N]) %*% Z %*% P[,,N])
  
  # return the error terms and the estimates
  list(a=a, P=P, at=at, Pt=Pt, v=v, G=G)
}

# forecast future states based on the Kalman filter estimates
kf_forecast <- function(kf_out, m) {
  # set up storage for projections
  a <- matrix(0, nrow=m+1, ncol=3)
  P <- array(NA, dim=c(3, 3, m+1))
  a[1,] <- kf_out$prev.mean
  P[,,1] <- kf_out$prev.var
  # generate projections
  for (i in 2:(m+1)) {
    a[i,] <- (kf_out$S %*% a[i-1,])
    P[,,i] <- kf_out$S %*% P[,,i-1] %*% t(kf_out$S) + (kf_out$R %*% kf_out$Q %*% t(kf_out$R))
  }
  # return mean and variance
  mn <- apply(a, 1, function(x) kf_out$Z %*% x)
  vr <- apply(P, 3, function(M) kf_out$Z %*% M %*% t(kf_out$Z))
  list(mean=mn, vr=vr)
}
