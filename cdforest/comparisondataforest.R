### Comparison Data Forest

# load packages

library(psych)
library(mvtnorm)
library(ineq)
library(purrr)
library(mlr) 

comparison.data.forest <- function(data, kmax, reps = 100, popsize = 10000, use = "pairwise.complete.obs", Cor.Type = "pearson"){
  
  # Generate comparison data
  
  comp.data <- matrix(nrow = kmax*reps, ncol = 183)
  
  for(i in 1:kmax){
    
    dt <- GenData(data, N.Factors = i, N = popsize, Cor.Type = Cor.Type, use = use)
    
    for(j in 1:reps){
      
      Sample <- dt[sample(1:popsize, size = nrow(data), replace = T),]
      
      # Calculate Features
      
      N <- nrow(Sample)
      p <- ncol(Sample)
      dat_cor <- cor(Sample, use = use)
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
      Gini <- ineq(dat_cor[lower.tri(dat_cor)], type = "Gini")
      Kolm <- ineq(dat_cor[lower.tri(dat_cor)], type = "Kolm")
      
      # reduced eigenvalues based on factor model
      
      fa_eigval <- fa(Sample, nfactors = 1, rotate="none", fm = "ml", warnings = FALSE)$values
      
      # empirical kaiser criterion
      
      ekc <- EKC(Sample)
      
      # setting missing eigenvalues to -1000
      
      eigval[(length(eigval)+1):80] <- -1000
      fa_eigval[(length(fa_eigval)+1):80] <- -1000
      names(eigval) <- paste("eigval", 1:80, sep = "")
      names(fa_eigval) <- paste("fa_eigval", 1:80, sep = "")
      
      # Fill in Comp.Data
      
      comp.data[reps*(i-1) + j,] <- unlist(cbind(data.frame(i, N, p, eiggreater1, releig1, releig2, releig3, eiggreater07,
                                                            sdeigval, var50, var75, onenorm, frobnorm, maxnorm, avgcor,
                                                            specnorm, smlcor, avgcom, det, KMO, Gini, Kolm, ekc), t(eigval),t(fa_eigval)))
      
    }
  }
  
  comp.data <- as.data.frame(comp.data)
  colnames(comp.data) <- c("k", "N", "p", "eiggreater1", "releig1", "releig2", "releig3", "eiggreater07",
                           "sdeigval", "var50", "var75", "onenorm", "frobnorm", "maxnorm", "avgcor", 
                           "specnorm", "smlcor", "avgcom", "det", "KMO", "Gini", "Kolm", "ekc",
                           paste("eigval", 1:80, sep = ""), paste("fa_eigval", 1:80, sep = ""))
  comp.data$k <- factor(comp.data$k)
  
  Task = makeClassifTask(data = comp.data, target = "k")
  forest = train("classif.ranger", Task)
  
  # Calculate Features for Empirical Data
  
  N <- nrow(data)
  p <- ncol(data)
  dat_cor <- cor(data, use = use)
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
  Gini <- ineq(dat_cor[lower.tri(dat_cor)], type = "Gini")
  Kolm <- ineq(dat_cor[lower.tri(dat_cor)], type = "Kolm")
  
  # reduced eigenvalues based on factor model
  
  fa_eigval <- fa(data, nfactors = 1, rotate="none", fm = "ml", warnings = FALSE)$values
  
  # empirical kaiser criterion
  
  ekc <- EKC(data)
  
  # setting missing eigenvalues to -1000
  
  eigval[(length(eigval)+1):80] <- -1000
  fa_eigval[(length(fa_eigval)+1):80] <- -1000
  names(eigval) <- paste("eigval", 1:80, sep = "")
  names(fa_eigval) <- paste("fa_eigval", 1:80, sep = "")
  
  # combination of features
  
  features <- cbind(data.frame(N, p, eiggreater1, releig1, releig2, releig3, eiggreater07, sdeigval,
                               var50, var75, onenorm,frobnorm, maxnorm, avgcor, specnorm, smlcor, avgcom,
                               det, KMO, Gini, Kolm, ekc), t(eigval), t(fa_eigval))
  
  predict(forest, newdata = features)$data$response
}