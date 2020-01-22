
########## Appends an object to a list #####################################
# list: is the list object, function appends new object at end of list     #
#  add: is the object to append to list                                    #
# name: is the name of this element object added                           #
############################################################################
appendlist <- function(list, add, name){
  length = length(list)+1                             # Length of the list, plus 1 for the new object added
  nameold = names(list)                               # Original names of list objects
  name = c(nameold, name)                             # Merging the old names with the new name of the new element
  list[[length]] <- add                               # Adding the new element to the list
  names(list) <- name                                 # Naming the new element
  return(list)                                        # Returning the newly appended list
}



########## Computes MVN Density values ####################################
#     x: is the values of x to compute the density at                     #
#    mu: is the vector of means of the random variables                   #
# sigma: is the Variance-Covariance matrix of the random variables        #
###########################################################################
dmvnorm = function(x, mu, sigma){
  #-    Input variables 
  #-
  #-    mu: vector of length p
  #- sigma: matrix of dim p x p
  #-     x: matrix of dim n x p OR vector of length p
  
  #-   Code variables
  #-
  #-   d: dimensions of x
  #-   p: dimensions of mean
  #- det: determinant of sigma matrix
  
  # Checking x is either a matrix, a vector, or a data.frame
  # Converts to matrix if data.frame
  if(!is.vector(x) && !is.matrix(x)){
    if(is.data.frame(x)){
      x = data.matrix(x)
    } else {
      stop("Error: x must be either a vector or a matrix")
    }
  }
  
  # Storing d for use later on
  if(is.vector(x)){
    d = length(x)
  } else {
    if(is.matrix(x)){
      d = ncol(x)
    }
  }
  
  # Defaulting to standardized means if none provided
  if(missing(mu)){
    mu = rep(0,d)
  }
  
  # Defaulting to standardized covariance if none provided
  if(missing(sigma)){
    sigma = diag(d)
  }
  
  # Checking means are vectors and taking dimensions for MVN
  # Also storing p for use later on
  if(!is.vector(mu)){
    stop("Error: mu needs to be a vector of length p")
  }else{
    p = length(mu)
  }
  
  # Checking sigma is a valid variance-covariance matrix
  # Also storing det for use later on
  if(!is.matrix(sigma) || !isSymmetric(sigma)){
    stop("Error: sigma must be a symmetric p by p positive-definite matrix")
  } else {
    if(ncol(sigma) != p){
      stop("Error: mu and sigma do not have conforming dimensions")
    } else {
      det = det(sigma)
      if(det<=0){
        stop("Error: sigma is not positive-definite")
      }
    }
  }
  
  # Checking if x is a vector or matrix of valid dimensions
  if(d != p){
    stop("Error: dimensions of x do not conform with mu and sigma")
  }
  
  # Returning probability of value if vector, or vector of probabilities if matrix
  # Utilizies solve() to determine inverse of sigma
  if(is.vector(x)){
    numer = exp(-0.5 * (x - mu) %*% solve(sigma) %*% (x-mu) )
    denom = sqrt( (2*pi)^p * det )
    return(numer/denom)
  } else {
    if(is.matrix(x)){
      # Trick to avoid loops, utilizes recycling and R's method of storing matrices
      v = t(x) - mu
      numer = exp(-0.5 * colSums( (solve(sigma) %*% v) * v))
      denom = sqrt( (2*pi)^p * det )
      return(numer/denom)
    } else{
      stop("Error: Something went very wrong")
    }
  }
}



########## Computes MVN Random values #####################################
#     n: is th number of pulls from the MVN distribution to pull          #
#    mu: is the vector of means of the random variables                   #
# sigma: is the Variance-Covariance matrix of the random variables        #
###########################################################################
rmvnorm = function(n, mu, sigma){
  #-    Input variables 
  #-
  #-     n: a numeric integer
  #-    mu: vector of length p
  #- sigma: matrix of dim p x p
  
  #-   Code variables
  #-
  #-   p.dim: dimensions of sigma (specifically columns)
  #-   d.dim: length of mu
  
  # Running checks on n
  if(length(n) != 1 || !is.numeric(n)){
    stop("Error: n must be a single number of pulls desired")
  } else {
    if(n%%1 != 0 || n <= 0){
      stop("Error: n must be a natural number")
    }
  }
  
  # Using standardized bivariate params if none given
  # Running standardized params of equiv dim if one is missing but other is provided
  if(missing(mu) || missing(sigma)){
    if(missing(mu) && missing(sigma)){
      mu = rep(0, 2)
      sigma = diag(2)
    }else{
      if(missing(mu) && !missing(sigma)){
        p.dim = ncol(sigma)
        mu = rep(0, p.dim)
      }else{
        d.dim = length(mu)
        sigma = diag(d.dim)
      }
    }
  }
  
  # Ensuring mu is a vector
  if(!is.vector(mu)){
    stop("Error: mu must be a vector")
  }
  
  # Ensuring sigma is a positive-definite symmetric matrix
  if(!is.matrix(sigma)){
    stop("Error: sigma must be a matrix")
  } else {
    if(!isSymmetric(sigma)){
      stop("Error: sigma must be symmetric")
    } else {
      if(det(sigma)<=0){
        stop("Error: sigma must be positive-definite")
      }
    }
  }
  
  # Taking dimensions of mu and sigma, if not done already
  if(!exists("p.dim")){
    p.dim = ncol(sigma)
  }
  if(!exists("d")){
    d.dim = length(mu)
  }
  
  # Checking mean dimensions match covariance matrix dimensions
  if(d.dim != p.dim){
    stop("Error: mu and sigma dimensions do not conform")
  }
  
  # Running affine transformation to generate values
  z <- matrix(rnorm(n*p.dim), nrow = n, ncol = p.dim)         # Initialize z ~ MVN(0,I)
  eig <- eigen(sigma, symmetric = TRUE)                       # Doing Eigendecomposition of Sigma
  U <- eig$vectors                                            # Eigenvectors
  V <- diag(sqrt(eig$values))                                 # Sqrt of Eigenvalues
  B <- U %*% V                                                # Sigma = B %*% t(B) = U %*% t(V) %*% U
  B <- t(B)                                                   # This transpose follows from transposing the theory; note Z  defined transposed to theory
  x <- z %*% B                                                # Applying Affine Transformation and returning mean zero MVN pulls
  x <- sweep(x, mu, MARGIN = 2, FUN = "+")                    # Adding means
  return(x)
}



########## Generate MVN Mixture data #######################################
#      n: is the integer number of pulls from the MVMN distribution        #
#    mus: is a list of vectors of means of the random variables            #
#   covs: is a list of Variance-Covariance matrices                        #
#  probs: is a vector of probabilities for each MVN                        #
# nclust: is the number of clusters intended                               #
############################################################################
rmvmixnorm = function(n, mus, covs, probs, nclust){
  #-    Input variables 
  #-
  #-      n: a natural number
  #-    mus: list of vectors each of length p
  #-   covs: list of matrices each of dim p x p
  #-  probs: a numeric integer
  #- nclust: a natural number greater than 2
  
  #-   Code variables
  #-
  #-    muslistlength: length of mus list
  #-  covslistlength: length of cos list
  #- probslistlength: length of probs list
  #-          mudims: length of mu vectors
  #-        covsdims: dimension of cov matrices
  
  # Running argument checks
  if(missing(n) ||
     missing(mus) ||
     missing(covs) ||
     missing(probs) ||
     missing(nclust)){
    stop("Error: Missing arguments")
  }
  # Checking n
  if(length(n) != 1 || !is.numeric(n)){
    stop("Error: n must be a single number of pulls desired")
  } else {
    if(n%%1 != 0 || n <= 0){
      stop("Error: n must be a natural number greater than 1")
    }
  }
  
  # Checking mus
  if(!is.list(mus)){
    stop("Error: mus must be passed as a list")
  } else {
    muslistlength = length(mus)
    if(muslistlength < 2){
      stop("Error: mus must be multiple vectors")
    } else {
      if(!all(sapply(mus, is.vector))){
        stop("Error: list of mus must be of vectors")
      } else {
        mudims = unique(sapply(mus, length))
        if(length(mudims) != 1){
          stop("Error: all mus in list must be vectors of same lengths")
        }
      }
    }
  }
  
  # Checking covs
  if(!is.list(covs)){
    stop("Error: covs must be a list")
  } else {
    covslistlength = length(covs)
    if(covslistlength < 2){
      stop("Error: covs must be multiple matrices")
    } else{ 
      if(!all(sapply(covs, is.matrix))){
        stop("Error: list of covs must be of matrices")
      } else {
        covdims = unique(as.vector(sapply(covs,dim)))
        if(length(covdims) != 1){
          stop("Error: all covariance matrices must have same dimensions AND be symmetric")
        } else {
          if(!all(sapply(covs,det) > 0)){
            stop("Error: all covs must be positive-definite")
          }
        }
      }
    }
  }
  
  # Checking probs
  if(!is.vector(probs)){
    stop("Error: probs must be a vector")
  } else {
    if(!all(probs > 0)){
      stop("Error: probs values must be positive")
    } else {
      if(!all(probs < 1)){
        stop("Error: probs values must be less than 1")
      } else {
        if(sum(probs) != 1){
          stop("Error: probs must sum to 1")
        } else{
          probslength = length(probs)
          if(probslength < 2){
            stop("Error: multiple probabilities must be given")
          }
        }
      }
    }
  }
  
  if(nclust != muslistlength){
    stop("Error: Number of clusters and number of mean vectors do not conform")
  } else {
    if(nclust != covslistlength){
      stop("Error: Number of clusters and number of covariance matrices do not match")
    } else {
      if(nclust != probslength){
        stop("Error: Number of clusters and number of probabilties do no match")
      }
    }
  }
  
  
  # Doing sampling from clusters
  coinflip = sample(c(1:nclust),
                    size = n,
                    prob = probs,
                    replace = TRUE)
  # Initialize output
  out = NULL
  # Computing pulls from Mixed MVN
  for(k in 1:nclust){
    noutcomes = sum(coinflip == k)
    if(noutcomes != 0){
      mean = mus[[k]]
      cov  = covs[[k]]
      out = rbind(out, rmvnorm(noutcomes, mu = mean, sigma = cov))
    }
  }
  # Shuffling output to avoid ordering created
  out = out[sample(nrow(out)), ]
  return(out)
}
















########## EM for gaussian mixture ##########
# 
#
###############################################################################################
  # Create starting values                                                                    
  # mustart = list(c(3, 60), c(3, 60.1))
  # covstart = list(cov(faithful), cov(faithful))
  # probs = c(.01, .99)
  # 
  # starts = list(mu=mustart, var=covstart, probs=probs)  # params is a list of mu var and probs
  # test = gaussmixEM(params=starts,
  # X=as.matrix(faithful),
  # clusters = 2, tol=1e-10,
  # maxits=1500)
  # http://georgepavlides.info/expectation-maximization-gaussian-mixtures-vectorized-matlab-octave-approach/
##############################################################################################
gaussmixEM = function(x, mus, covs, probs, clusters = 2, tol = 1e-10, iter = 1e3){
  
  n <- nrow(x)                                         # Counting number of points
  d <- ncol(x)                                         # Counting dimensions
  
  # # Checking valid parameters
  # if(!sum(probs) == 1 ||
  #    !all(probs > 0) ||
  #    !all(lapply(covs,det)>=0) ||
  #    !all(sapply(covs,isSymmetric)) ||
  #    !length(mus) == length(covs) ||
  #    !length(mus) == clusters ||
  #    !clusters >= 2 ||
  #    !d == length(mus[[1]]) ||
  #    !is.list(covs) ||
  #    !is.list(mus) ){
  #   stop("Argument error")
  # }
  # 
  # # Checking included parameters
  # if(missing(x) ||
  #    missing(mus) ||
  #    missing(covs) ||
  #    missing(probs)){
  #   stop("Missing arguments")
  # }
  
  # Initialization
  outlist <- as.list(0)                                # Initializing output
  converged <- FALSE                                   # Initializing convergence check
  count <- 0                                           # Initializing iteration counts
  #loglik <- 0                                         # Initializing loglikelihood value
  res <- matrix(0, nrow = n, ncol = clusters)          # Initializing responsibilities; notice transposed from expected
  
  # Initial loglikelihood
  loglik <- matrix(0, nrow = n, ncol = clusters)
  for(k in 1:clusters){
    loglik[,k] <- t(probs[k]) %*% dmvnorm(x, mus[[k]],covs[[k]])
  }
  loglik <- log(rowSums(loglik))
  loglik <- sum(loglik)
  
  # Looping EM Steps until convergence
  while(count < iter || !converged){
    probs.old <- probs
    loglik.old <- loglik
    res.old <- res
    
    ## E Step
    for(k in 1:clusters){
      res[ , k] <- probs[k] * dmvnorm(x,
                                     mus[[k]],
                                     covs[[k]])
    }
    res = res / rowSums(res)                          # Responsibilities computed
    res.sum = colSums(res)                            # Sum of cluster responsibilities; used for parameters
    ## M Step
    probs <- colMeans(res)                            # New probabilities computed
    for(k in 1:clusters){
      mus[[k]] <- colSums(res[,k] * x) / res.sum[k]            # New means computed; note the numerator is transposed from expected
    }
    #return(list(mu = mus, probs = probs, res = res, res.sum = res.sum)) # used for checking
    for(k in 1:clusters){
      cov <- matrix(0, nrow = ncol(x), ncol = ncol(x))
      mu <- as.vector(drop(mus[[k]]))
      for(i in 1:n){
        x <- as.matrix(x[i,])
        covs[[k]] <- colSums(res[,k] * (x - mu) %*% t((x - mu))) / res.sum[k]
        # # y <- as.matrix(x[i,])
        # # cov <- cov + res[i,k] * y%*%t(y)
        # y <- as.matrix(x[i, ])
        # cov = cov + res[i,k] * (y - drop(mus[[k]])) %*% t(y - drop(mus[[k]]))
      }
      # covs[[k]] <- cov[[k]] / res.sum
      # cov[[k]] <- cov / res.sum[k] - mus[[k]] %*% t(mus[[k]])
    }
    #return(list(mu = mus, probs = probs, res = res, res.sum = res.sum)) # used for checking
    
    # loglik <- matrix(0, nrow = n, ncol = clusters)
    # for(k in 1:clusters){
    #   for(i in 1:n){
    #     loglik[,k] <- t(probs[k]) %*% dmvnorm(x, mus[[k]],covs[[k]])
    #   }
    # }
    # loglik <- log(rowSums(loglik))
    # loglik <- sum(loglik)
    # converged <- (loglik - loglik.old) <= tol
    # count <- count + 1
    # if(converged || count == iter){
    #   return(list(mus = mus, covs = covs, probs = probs, iterations = count))
    # }
  }
}



########## Plot EM ##########
plot.gaussmixEM <- function(data, object, clusters, which.iter = object[[1]],
                            xlabels = NULL, ylabels = NULL,
                            xlimits = c(min(data[,1]),max(data[,1])),
                            ylimits = c(min(data[,2]),max(data[,2]))){
  expmax <- object[[which.iter]]
  plot(data, type = "n",
       xlab = xlabels, ylab = ylabels,
       xlim = xlimits, ylim = ylimits)
  for(j in 1:clusters){
    points(data[expmax$cluster == j, ], pch = j)
    x.points <- seq(0.5*expmax$mu[j,1], 2*expmax$mu[j,1], length.out = 100)
    y.points <- seq(0.5*expmax$mu[j,2], 2*expmax$mu[j,2], length.out = 100)
    z <- matrix(0, nrow = 100, ncol = 100)
    mu <- expmax$mu[j,]
    sigma <- expmax$var[[j]]
    for (x in 1:100) {
      for (y in 1:100) {
        z[x,y] <-dmvnorm(c(x.points[x],  y.points[y]), mean=mu, sigma=sigma)
      }
    }
    contour(x.points,y.points,z, add = TRUE)
  }
}