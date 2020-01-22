## These scripts are badly organized.

EMmmvn = function(x, mus, covs, probs, nclust = 2, tol = 1e-10, iter = 1e3){
  #- Input variables
  #-      x: A matrix where each column is a variable, and each row is an observation
  #-   ***NOTICE*** (Order matters for mus, covs, and probs; i.e. mus[2,] and covs[,,2] and prob[2] are params for cluster 2)
  #-    mus: A matrix where each row is a mean vector of a cluster  
  #-   covs: A  3-d array of matrices of Variance-Covariance matrices of a cluster 
  #-  probs: A vector of probabilities for each cluster
  #- nclust: The number of inteded clusters (note that this could be unnecessary, but is used for several sanity checks)
  #-    tol: A small number used to check for convergence from one iteration to the next
  #-   iter: The maximum number of possible iterations, when count reaches this number the algorithm stops
  
  #- Code Variables
  #-          N: Number of Observations/"Points" available
  #-          D: Dimensions of Mixed-MVN / Number of Variables
  #-  converged: A boole to check whether convergence requirements have been met (initial:FALSE)
  #-        res: A matrix of responsibilities for each Observation/"Point"
  #-      count: A counter for the number of iterations made 
  
  # Example operation
  #       x = data.matrix(faithful)
  #     mus = rbind(c(3,59), c(3,60))
  #    covs = array(cov(faithful), dim = c(2,2,2))
  #   probs = c(0.01, 0.99)
  #  nclust = 2
  #     tol = 1e-10
  #    iter = 1e3
  #  output = EMmmvn(x, mus, covs, probs, nclust, tol, iter)
  
  N <- nrow(x)
  D <- ncol(x)
  
  # Initilization
  outlist <- list(0)
  converged <- FALSE
  count <- 0
  res <- matrix(0, nrow = nclust, ncol = N)
  loglik <- matrix(0,nrow = nclust, ncol = nrow(x))
  for(k in 1:nclust){
    loglik[k,] = probs[k]*EMmmvnDMVNORM(x, mus[k,], covs[,,k], loglike = FALSE)
  }
  loglik = log(colSums(loglik))
  loglik = sum(loglik)
  # Doing updates
  while( !converged ){
    # Storing previous values for comparison
    mus.old <- mus
    covs.old <- covs
    probs.old <- probs
    loglik.old <- loglik
    res.old <- res
    
    # --- E Step ---
    for(k in 1:nclust){
      res[k,] = probs[k] * EMmmvnDMVNORM(x, mus[k, ], covs[,,k], loglike = FALSE)
    }
    res = sweep(res, 2, colSums(res), FUN = "/")
    res.sum = rowSums(res)
    # --- M step ---
    mus = res %*% x / res.sum
    probs = rowMeans(res)
    loglik <- matrix(0,nrow = nclust, ncol = nrow(x)) # DELETE
    for(k in 1:nclust){                                               # Covariance doesn't simplify?
      cov = matrix(0, D, D)
      for(i in 1:N){
        cov = cov + res[k, i] * tcrossprod(x[i,] - mus[k, ])
      }
      cov = cov / sum(res[k, ])
      covs[ , , k] = cov
      # Also handling loglikelihood in this loop
      loglik[k,] = probs[k]*EMmmvnDMVNORM(x, mus[k,], covs[,,k], loglike = FALSE)
    }
    loglik = log(colSums(loglik))
    loglik = sum(loglik)
    
    count = count + 1
    
    check.mus = min(abs(mus-mus.old)) <= tol
    check.covs = min(abs(covs-covs.old)) <= tol
    check.probs = min(abs(probs-probs.old)) <= tol
    check.loglike = min(abs(loglik-loglik.old)) <= tol
    
    converged = check.mus || check.covs || check.probs || check.loglike || count == iter
  }
  # Forming output
  outlist[[1]] = count
  outlist[[2]] = mus
  outlist[[3]] = covs
  outlist[[4]] = t(res)
  outlist[[5]] = loglik
  # This creates cluster memberships by assigning each row to the cluster it has the largest responsibility in
  outlist[[6]] = rowSums(sweep(t(apply(t(res), 1, function(x) as.numeric(x == max(x)))),c(1:nclust), MARGIN = 2,FUN = "*"))
  names(outlist) <- c("Iterations",
                      "Means",
                      "Covariances",
                      "Responsibilities",
                      "LogLikelihood",
                      "Cluster")
  return(outlist)
}


EMmmvnDMVNORM = function(x, mu, sigma, loglike = TRUE){
  # This is a modification of the dmvnorm function to be more efficient for use in the E-M algorithm
  p = length(mu)
  det = determinant(sigma,logarithm = TRUE)
  det = prod(unlist(det))
  inv.sigma = qr.solve(sigma) #chol2inv(chol(sigma, pivot = FALSE))#solve(sigma)  
  # qr.solve() is most reliable when near degenerate
  if(is.vector(x)){
    retval = -0.5 * (x-mu) %*% inv.sigma %*% (x - mu) -0.5 * det
  } else {
    if(is.matrix(x)){
      # Trick to avoid loops, utilizes recycling and R's method of storing matrices
      v = t(x) - mu
      retval = -0.5 * colSums( (inv.sigma %*% v) * v) - 0.5 * det
    } else{
      stop("Error: Something went very wrong") # Shouldn't flag unless an array is passed
    }
    if(loglike){
      return(retval)
    } else {
      expretval = exp(retval - p/2 * log(2*pi))
      return(expretval)
    }
  }
}

ExpectationMaximization = function(x, nclust, tol = 1e-10, iter = 1e3, nsearch = 100){
  ll = 0
  EM = list(0)
  for(i in 1:nsearch){
    tryCatch({
      mus = x[sample(nrow(x),nclust),]
      covs = array(rep(cov(x),nclust), dim = c(ncol(x),ncol(x),nclust))
      probs = rep(1/nclust,nclust)
      EM.result = EMmmvn(x = x , mus = mus,
                         covs = covs, probs = probs,
                         nclust = nclust, tol = tol, iter = iter)
      ll[i] = EM.result$LogLikelihood
      EM[[i]] = EM.result
    }, error = function(e){
      ll[i] = NA
      EM[[i]] = NA
    })
  }
  ll[ll == Inf] <- NA
  ll[ll == -Inf] <- NA
  result = EM[[which.max(ll)]]
  return(result)
}

#### Example 1 iris data
#### HIGHLY recommended to jitter data, otherwise degenerate covariances come up
# x = data.frame(lapply(iris[,-5], jitter))
# x = data.matrix(x)
# nclust = 3
# test = ExpectationMaximization(x, nclust, nsearch = 100)
# Cluster = as.factor(test$Cluster)
# levels(Cluster) = paste0("Cluster ",1:nclust)
# x = as.data.frame(x)
# x.em = cbind(x, Species = iris$Species, Cluster = Cluster)
# plot(Sepal.Length ~ Petal.Length,
#      pch = c(16,17,18)[Species],
#      col = (rep(1:nclust)+1)[Cluster],
#      data = x.em)
# table = table(x.em$Species,x.em$Cluster)
# print(table)
# print(test$LogLikelihood)
#pairs(x.em[,1:4], pch = c(16,17,18)[x.em$Species],col = (rep(1:nclust)+1)[x.em$Cluster])


# #### Example 2 faithful data
# #### Pivot: FALSE for best results
# #### qr.solve: Works
# x = data.matrix(faithful)
# nclust = 2
# test = ExpectationMaximization(x, nclust, nsearch = 10)
# Cluster = as.factor(test$Cluster)
# levels(Cluster) = paste0("Cluster", 1:nclust)
# x = as.data.frame(x)
# x.em = cbind(x, Cluster = Cluster)
# plot(waiting ~ eruptions,
#      col = (rep(1:nclust)+1)[Cluster],
#      data = x.em)

# x = jitter(data.matrix(iris[,-5]))
# Factor = iris$Species
# nclust = 3
# test = ExpectationMaximization(x, nclust)
# x = as.data.frame(x)
# Cluster = as.factor(test$Cluster)
# levels(Cluster) = paste("Cluster",1:nclust)
# x.em = cbind(x, Factor, Cluster)
# plot(data = x.em, Sepal.Length~Petal.Length, col = (seq(1:nclust)+1)[x.em$Cluster], pch = (seq(1:nlevels(Factor))+16)[x.em$Factor])
# Cluster = as.factor(test$Cluster)
# levels(Cluster) = paste0("Cluster ",1:nclust)
# print(table(iris$Species, Cluster))

# plot_ly(data = x.em, y = ~Sepal.Length, x = ~Petal.Length, mode = "markers", type = "scatter", color = ~Cluster, symbol = ~Factor)



