#' Run Step 1 of Multivariate Outcome Regression method proposed by Sofer et al (2014)
#'
#' @param Xaug Design matrix of original data column augmented by knockoffs
#' @param Y vectorized response vector
#' @param n number of samples
#' @param K number of responses
#' @return the omega estimated through step 1
#' @export Step1
Step1 = function(Xmat, Y, n, K, method = "lasso"){
  Omega0 = diag(K)

  if(method == "scad"){
    sel = SCADSelection(Xmat, Y, n, Omega0)
  }
  if(method == "lasso"){
    sel = LassoSelection(Xmat, Y, n, Omega0)
  }
  if(method == "alasso"){
    sel = adaLassoSelection(Xmat, Y, n, Omega0)
  }
  optBeta = sel$optimalBeta
  resids = Y-Xmat%*%t(t(optBeta))
  residMatrix = matrix(0, nrow = n, ncol = K)
  for(i in 1:K){
    marg = seq(i,length(resids), by = K)
    residMatrix[,i] = marg
  }

  glassoStep1 = CVglasso(X=residMatrix)
  step1Omega = glassoStep1$Omega
  return(step1Omega)
}

#' Run Step 2 of Multivariate Outcome Regression method proposed by Sofer et al (2014)
#'
#' @param Xaug Design matrix of original data column augmented by knockoffs
#' @param Y vectorized response vector
#' @param step1Omega the omega estimated through step 1
#' @param n the number of samples
#' @param penalty the penalty function of interest
#' @return The result after selection
#' @export Step2
Step2 = function(Xaug, Y, step1Omega, n, method = "scad"){
  if(method == "scad"){
    result = SCADSelection(Xaug, Y, n, step1Omega)
  }
  if(method == "lasso"){
    result = LassoSelection(Xaug, Y, n, step1Omega)
  }
  if(method == "alasso"){
    result = adaLassoSelection(Xaug, Y, n, step1Omega)
  }
  return(result)
}

#' Run the selection method through Lasso penalty (L1)
#'
#' @param Xaug Design matrix of original data column augmented by knockoffs
#' @param Y vectorized response vector
#' @param n the number of samples
#' @param Omega0 the estimated omega matrix
#' @return The result of the lasso fit
LassoSelection = function(Xaug, Y, n, Omega0 = NULL){
  Xaug = scale(Xaug)
  lassoFit = cv.glmnet(Xaug, Y, intercept = F)
  optimalBeta = coef(lassoFit, s = "lambda.min")[-1]
  res = list(fit = lassoFit, optimalBeta = optimalBeta)
  return(res)
}

#' Run the selection method employing the adaptive lasso penalty (Zou, 2006)
#'
#' @param Xaug Design matrix of original data column augmented by knockoffs
#' @param Y vectorized response vector
#' @param n the number of samples
#' @param Omega0 the estimated omega matrix
#' @return The result of adaptive lasso method
#' @export adaLassoSelection
adaLassoSelection = function(Xaug, Y, n, Omega0 = NULL){
  Xaug = scale(Xaug)
  ridge1 = cv.glmnet(x = Xaug, y = Y, alpha = 0)
  ridgeCoef = coef(ridge1, s = ridge1$lambda.min)[-1]

  alasso1 = cv.glmnet(Xaug, Y, alpha = 1, penalty.factor = 1/abs(ridgeCoef))
  optimalBeta = coef(alasso1, s = "lambda.min")[-1]
  res = list(fit = alasso1, optimalBeta = optimalBeta)
  return(res)
}

#' Run the selection method employing the SCAD penalty (Fan, 2001)
#'
#' @param Xaug Design matrix of original data column augmented by knockoffs
#' @param Y vectorized response vector
#' @param n the number of samples
#' @param Omega0 the estimated omega matrix
#' @return The result from the SCAD fit
#' @export SCADSelection
SCADSelection = function(Xaug, Y, n, Omega0 = NULL){
  # Run SCAD
  scad.l2 = cv.ncvreg(Xaug,  Y, penalty = "SCAD")
  optimalBeta = coef(scad.l2, s = "lambda.min")[-1]
  scadFit = scad.l2
  res = list(optimalBeta = optimalBeta, scadFit = scadFit)
  return(res)
}

# Single Knockoff Setting
#' @export runSingleKF
runSingleKF = function(X,Y,q, penalty){
  if(penalty == "lasso"){
    print("Running Lasso Stat")
    penaltyStat = stat.glmnet_coefdiff
  }
  if(penalty == "scad"){
    print("Running SCAD Stat")
    penaltyStat = stat.SCAD_coefdiff
  }
  if(penalty == "alasso"){
    print("Running A.Lasso Stat")
    penaltyStat = stat.ada_coefdiff
  }


  kfp = knockoff.filter(X, Y, statistic = penaltyStat, fdr = q, offset = 1)$selected
  kf = knockoff.filter(X, Y, statistic = penaltyStat, fdr = q, offset = 0)$selected
  res = list(S = kf, Sp = kfp)
  return(res)

}


#' @export stat.SCAD_coefdiff
stat.SCAD_coefdiff <- function(X, X_k, y, family='gaussian', cores=2, ...) {
  # Compute statistics
  Z = cv_coeffs_SCAD(cbind(X, X_k), y, family=family, parallel=parallel, ...)
  p = ncol(X)
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[orig+p])

}
#' @export cv_coefffs_SCAD
cv_coeffs_SCAD <- function(X, y, nlambda=500, intercept=T, parallel=T, ...) {
  n = nrow(X); p = ncol(X)
  cv.SCAD = cv.ncvreg(X, y, penalty = "SCAD")
  coef(cv.SCAD, s = "lambda.min")[2:(p+1)]
}

#' @export stat.ada_coefdiff
stat.ada_coefdiff <- function(X, X_k, y, family='gaussian', cores=2, ...) {
  # Compute statistics
  Z = cv_coeffs_ada(cbind(X, X_k), y, family=family, parallel=parallel, ...)
  p = ncol(X)
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[orig+p])

}
#' @export cv_coeffs_ada
cv_coeffs_ada <- function(X, y, nlambda=500, intercept=T, parallel=T, ...) {
  n = nrow(X); p = ncol(X)

  ridge1 = cv.glmnet(x = X, y = y, alpha = 0)
  ridgeCoef = coef(ridge1, s = ridge1$lambda.min)[-1]

  alasso1 = cv.glmnet(X, y, alpha = 1, penalty.factor = 1/abs(ridgeCoef))
  coef(alasso1, s = "lambda.min")[2:(p+1)]
}
