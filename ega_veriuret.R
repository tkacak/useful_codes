# Generating the data:

library(tictoc)

bigmat <- vector("list", 500)
mat <- vector("list")

for(i in 1:nrow(levels.all)){
  mat[[i]] <- i
}


for(j in 1:500){
  bigmat[[j]] <- mat
}

temp <- NULL

library(parallel)
## Step 1: Create a cluster of child processes 
cl <- makeCluster( 5 )

## Step 2: Load the necessary R package(s)
## N.B. length(cl) is the number of child processes
##      in the cluster 
par.setup <- parLapply( cl, 1:length(cl),
                        function(xx) {
                          require(Matrix)
                          require(EGAnet)
                          require(qgraph)
                          require(NetworkToolbox)
                          require(mvtnorm)
                        })

## Step 3: Distribute the necessary R objects 
clusterExport( cl, c('bigmat', 'mat', "data.gen", "EGA", "levels.all"))

library(tictoc)
## Step 4: Do the computation
tic()
for(i in 1:500){
  bigmat[[i]] <- parLapply(cl, mat,
                           function(x) {
                             repeat{
                               temp <- data.gen(nvar=levels.all[x,2], nfac=levels.all[x,3], n = levels.all[x,1], load = levels.all[x,4], corf = levels.all[x,5])
                               ega.res <- is.na(tryCatch(EGA(temp, plot = FALSE)$n.dim,
                                                         warn=FALSE, error=function(e) return(NA)))
                               if(ega.res == FALSE){
                                 break
                               }
                             }
                             return(temp)
                           })
}



toc()

stopCluster(cl)
