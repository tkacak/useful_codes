load_mat <- function(k, p, vpf, pf, sf,vpf_unique = TRUE){
  
  switch(pf,
         "low" = {
           pf_low = 0.35 
           pf_upper = 0.5
         },
         "medium" = {
           pf_low = 0.5
           pf_upper = 0.65
         },
         "high" = {
           pf_low = 0.65
           pf_upper = 0.8
         })
  switch(sf,
         "none" = {
           sf_low = 0 
           sf_upper = 0
         },
         "low" = {
           sf_low = 0
           sf_upper = 0.1
         },
         "medium" = {
           sf_low = 0.1
           sf_upper = 0.2
         })
  
  x <- runif(p, pf_low, pf_upper) 
  y <- runif(p*(k-1) , sf_low, sf_upper)
  
  i <- 1:(p)
  if(vpf_unique){
    j <- rep(1:k, each=vpf) 
  }
  else{
    j <- rep(1:k, times=vpf) 
  }
  
  L <- matrix(NA, p, k) 
  L[cbind(i, j)] <- x 
  L[is.na(L)] <- y
  L
}

#Generate ordinal datasets with normal and different two levels skewed distributions
ordinal_simulation <- function( N,vpf, k, rho, pf, sf, distr, cat){
  # simulate data
  # simulate loading pattern
  p = vpf*k
  L <- load_mat(k, p, vpf, pf, sf)
  fcor <- matrix(rho, k,k) + diag(1-rho,k,k)
  Sigma <- L%*%fcor%*%t(L) + diag(diag(diag(p)-L%*%fcor%*%t(L)))
  
  # simulate data with specific distribution depending on simulation condition
  
  dat <- switch(distr,
                "normal" = {
                  dat <- data.frame(mvtnorm::rmvnorm(N, mean = rep(0,p), sigma = Sigma))
                  dat <- data.frame(
                    lapply(dat,
                           FUN = function(x){cut(x, c(min(x)-0.01,
                                                      qnorm((1:(cat-1))/cat),
                                                      max(x) + 0.01), labels = FALSE,
                                                 right =FALSE)}))
                  dat
                },
                "moderately" = {
                  
                  dat <- data.frame(mvtnorm::rmvnorm(N, mean = rep(0,p), sigma = Sigma))
                  dat <- data.frame(
                    lapply(dat,
                           FUN = function(x){cut(x, c(min(x)-0.01,
                                                      sort(sample(c(-1,1),1)*qnorm((1:(cat-1))/(1.5*cat))),
                                                      max(x) + 0.01), labels = FALSE,
                                                 right =FALSE)}))
                  dat
                },
                "severely" = {
                  
                  dat <- data.frame(mvtnorm::rmvnorm(N, mean = rep(0,p), sigma = Sigma))
                  dat <- data.frame(
                    lapply(dat,
                           FUN = function(x){cut(x, c(min(x)-0.01,
                                                      sort(sample(c(-1,1),1)*qnorm((1:(cat-1))/(2.5*cat))),
                                                      max(x) + 0.01), labels = FALSE,
                                                 right =FALSE)}))
                  dat}
  )#go switch function
  list(dat = dat, L = L, p = p, k = k, N = N, distr = distr)
}