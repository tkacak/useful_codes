# This script contains the function that predicts number of factors
# based on 184 features that are calculated for the empirical data set
# that has to be specified for the newdata argument
# the trained model is passed to the function via the mod argument

# the argument "use" can be used just like the respective argument in the cor function


factor.forest <- function(newdata, mod, use = "pairwise.complete.obs"){
  # calculate all necessary features:
  
  N <- nrow(newdata)
  p <- ncol(newdata)
  dat_cor <- cor(newdata, use = use)#EGAnet::polychoric.matrix(veri_kategorik[[i]])
  eigval <- eigen(dat_cor)$values
  vareig <- cumsum(eigval)/p
  
  # eigenvalue features
  
  eiggreater1 <- sum(eigval > 1)  
  releig1 <- eigval[1]/p
  releig2 <- sum(eigval[1:2])/p
  releig3 <- sum(eigval[1:3])/p
  eiggreater07 <- sum(eigval > 0.7)
  sdeigval <- sd(eigval)
  var50 <- min(which(vareig > 0.50))
  var75 <- min(which(vareig > 0.75))
  
  # matrix norm features
  
  onenorm <- norm(dat_cor,"O")
  frobnorm <- norm(dat_cor,"F")
  maxnorm <- norm(dat_cor-diag(p),"M")
  avgcor <- sum(abs(dat_cor-diag(p)))/(p*(p-1))
  specnorm <- sqrt(eigen(t(dat_cor)%*%dat_cor)$values[1])
  
  smlcor <- sum(dat_cor <= 0.1)
  avgcom <- mean(smc(dat_cor))
  det <- det(dat_cor)
  
  KMO <- KMO(dat_cor)$MSA
  Gini <- ineq(lower.tri(dat_cor), type = "Gini")
  Kolm <- ineq(lower.tri(dat_cor), type = "Kolm")
  
  # parallel analysis
  
  quiet(pa <- fa.parallel(newdata,fa="fa",plot = FALSE))
  pa_solution <- pa$nfact
  fa_eigval <- pa$fa.values
  
  
  # empirical kaiser criterion
  
  lref <- rep(0,p)
  for (i in 1:p) {
    lref[i] <- max(((1 + sqrt(p/N))^2) * (p-sum(lref))/(p-i+1),1)
  }
  ekc <- which(eigval<=lref)[1]-1
  
  # setting missing eigenvalues to -1000
  
  eigval[(length(eigval)+1):80] <- -1000
  fa_eigval[(length(fa_eigval)+1):80] <- -1000
  names(eigval) <- paste("eigval", 1:80, sep = "")
  names(fa_eigval) <- paste("fa_eigval", 1:80, sep = "")
  
  
  cd <- EFA.Comp.Data(Data = newdata, F.Max = 8, use = use)
  
  # combination of features
  
  features <- cbind(data.frame(N,p,eiggreater1,releig1,releig2,releig3,eiggreater07,sdeigval,var50,var75,onenorm,frobnorm,
                               maxnorm, avgcor, specnorm, smlcor, avgcom,det, KMO, Gini, Kolm, pa_solution, ekc, cd), t(eigval), t(fa_eigval))
  
  out <- predict(mod, newdata = features)$data
  
  cat("\n\n", "Probabilities for different factor solutions:", "\n\n")
  print(round(out[-9], digits = 5))
  cat("\n\n", paste("Suggested number of factors:", out$response, sep =" "), "\n\n")
  return(data.frame(out$response))
}
