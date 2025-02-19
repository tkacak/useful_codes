#Generate continuous datasets
continuous_simulation <- function(N,vpf, k, rho, pf, sf){
  # simulate data mvnormal continuous data
  
  p = vpf*k
  
  L <- load_mat(k, p, vpf, pf, sf)
  fcor <- matrix(rho, k,k) + diag(1-rho,k,k)
  Sigma <- L%*%fcor%*%t(L) + diag(diag(diag(p)-L%*%fcor%*%t(L)))
  
  dat <- data.frame(mvtnorm::rmvnorm(N, mean = rep(3,p), sigma = Sigma))
  
  list(dat = dat, L = L, p = p, k = k, N = N)
}