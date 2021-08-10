### Simulation Functions
#' @export generateData
generateData = function(n, p, nGroup, errorRho, Bmatrix, sigma, type){
  pstar = p+1
  packageInitialization()
  means = rep(1, pstar)

  print("Generating Cov")

  # compound symmetry
  if(type == "CS"){
    errorCovariance = matrix(errorRho, nGroup, nGroup)
    diag(errorCovariance) = rep(1, nGroup)
  }

  if(type  == "AR"){
    errorCovariance = matrix(0,nGroup, nGroup)
    for(i in 1:nGroup){
      for(j in 1:nGroup){
        errorCovariance[i,j] = errorRho^abs(i-j)
      }
    }
  }
  print("Genreating E")
  E = matrix(0, nrow = n, ncol = nGroup)
  for(i in 1:n){
    ei = rnorm(nGroup, mean = rep(0, nGroup), errorCovariance)
    E[i,] = ei
  }

  print("Generating W")
  # Generate W
  W = mvlognormal(n, means, Sigma = diag(sigma), R = sigma)

  for(i in 1:nrow(W)){
    for(j in 1:ncol(W)){
      if(W[i,j]==0) {W[i,j]=0.5}
    }
  }
  print("Generating Z")
  # Compute Z
  Z = acomp(W)
  ref = Z[,pstar]

  print("Generating X")
  # log transformation
  X = matrix(0, nrow = n, ncol = p)
  for(t in 1:ncol(X)){
    X[,t] = log(Z[,t]/ref)
  }
  betaVec = as.vector(Bmatrix)

  Ymat = X%*%Bmatrix + E
  Y = as.vector(Ymat)

  res = list(Y=Y, Ymat = Ymat, X = X, B = Bmatrix, betaVec = betaVec, W = W)
  return(res)
}

#' @export computeFDR
computeFDR = function(X, truth){
  d = max(length(X),1)
  FP = length(X) - length(which(X%in%truth))
  return(FP/d)
}

#' @export computePower
computePower = function(X,truth){
  P= length(truth)
  TP = length(which(X%in%truth))
  return(TP/P)
}
