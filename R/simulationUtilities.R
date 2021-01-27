### Simulation Functions
generateData = function(n,p,q, nGroup, s1, s2,errorRho, Bmatrix,sigma){
  pstar = p+1
  loadPackages()
  means = rep(1, pstar)

  print("Generating Cov")
  # error covariance
  errorCovariance = matrix(0,nGroup, nGroup)
  for(i in 1:nGroup){
    for(j in 1:nGroup){
      errorCovariance[i,j] = errorRho^abs(i-j)
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

  print("Generating KF")
  # Generate Knockoffs
  Xkf = create.second_order(X)

  Im = diag(nGroup)

  for(i in 1:n){
    Xi = kronecker(t(X[i,]), Im)
    Xki = kronecker(t(Xkf[i,]), Im)
    Ei  = E[i,]
    tXi = t(Xi)
    tXki = t(Xki)
    tEi = t(E[i,])
    if(i == 1){
      Xmat = tXi
      Xkmat = tXki
      Emat = tEi
    }
    else{
      Xmat = cbind(Xmat, tXi)
      Xkmat = cbind(Xkmat, tXki)
      Emat = cbind(Emat, tEi)
    }
  }
  Xmat = t(Xmat)
  Xkmat = t(Xkmat)
  Emat = t(Emat)
  Xaug = cbind(Xmat, Xkmat)
  betaVec = as.vector(t(Bmatrix))
  Y = Xmat%*%betaVec + Emat

  res = list(Y=Y, X = X, B = Bmatrix, betaVec = betaVec, Xaug = Xaug, W = W)
  return(res)
}

computeFDR = function(X, truth){
  d = max(length(X),1)
  FP = length(X) - length(which(X%in%truth))
  return(FP/d)
}

computePower = function(X,truth){
  P= length(truth)
  TP = length(which(X%in%truth))
  return(TP/P)
}


#### Marginal SCAD
stat.SCAD_coefdiff <- function(X, X_k, y, family='gaussian', cores=2, ...) {
  # Compute statistics
  Z = cv_coeffs_SCAD(cbind(X, X_k), y, family=family, parallel=parallel, ...)
  p = ncol(X)
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[orig+p])

}
cv_coeffs_SCAD <- function(X, y, nlambda=500, intercept=T, parallel=T, ...) {
  n = nrow(X); p = ncol(X)
  cv.SCAD = cv.ncvreg(X, y, penalty = "SCAD")
  coef(cv.SCAD, s = "lambda.min")[2:(p+1)]
}


#### Marginal Adaptive
stat.ada_coefdiff <- function(X, X_k, y, family='gaussian', cores=2, ...) {
  # Compute statistics
  Z = cv_coeffs_ada(cbind(X, X_k), y, family=family, parallel=parallel, ...)
  p = ncol(X)
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[orig+p])

}
cv_coeffs_ada <- function(X, y, nlambda=500, intercept=T, parallel=T, ...) {
  n = nrow(X); p = ncol(X)

  ridge1 = cv.glmnet(x = X, y = y, alpha = 0)
  ridgeCoef = coef(ridge1, s = ridge1$lambda.min)[-1]

  alasso1 = cv.glmnet(X, y, alpha = 1, penalty.factor = 1/abs(ridgeCoef))
  coef(alasso1, s = "lambda.min")[2:(p+1)]
}
