#' Run Step 1 of Multivariate Outcome Regression method proposed by Sofer et al (2014)
#'
#' @param Xaug Design matrix of original data column augmented by knockoffs
#' @param Y vectorized response vector
#' @param n number of samples
#' @param K number of responses
#' @return the omega estimated through step 1
Step1 = function(Xaug, Y, n, K, method = "S"){
  Omega0 = diag(K)

  if(method == "S"){
    sel = SCADSelection(Xaug, Y, n, Omega0)
  }
  if(method == "L"){
    sel = LassoSelection(Xaug, Y, n, Omega0)
  }
  if(method == "A"){
    sel = adaLassoSelection(Xaug, Y, n, Omega0)
  }
  optBeta = sel$optimalBeta
  resids = Y-Xaug%*%t(t(optBeta))
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
Step2 = function(Xaug, Y, step1Omega, n, penalty = "scad"){
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
  for(i in 1:n){
    if(i == 1){
      blockOmega0 = Omega0
    }
    else{
      blockOmega0 = bdiag(blockOmega0, Omega0)
    }
  }
  sqrtOmega0 = sqrtm(blockOmega0)
  whiteY = sqrtOmega0%*%Y
  whiteXaug = scale(sqrtOmega0%*%Xaug)

  lassoFit = cv.glmnet(whiteXaug, whiteY, intercept = F)
  optimalBeta = coef(lassoFit, s = "lambda.min")[-1]
  res = list(fit = lassoFit, optimalBeta = optimalBeta)
  return(res)
}

#' Run the selection method through adaptive lasso penalty (Zou, 2006)
#'
#' @param Xaug Design matrix of original data column augmented by knockoffs
#' @param Y vectorized response vector
#' @param n the number of samples
#' @param Omega0 the estimated omega matrix
#' @return The result of adaptive lasso method
adaLassoSelection = function(Xaug, Y, n, Omega0 = NULL){
  for(i in 1:n){
    if(i == 1){
      blockOmega0 = Omega0
    }
    else{
      blockOmega0 = bdiag(blockOmega0, Omega0)
    }
  }
  sqrtOmega0 = sqrtm(blockOmega0)
  whiteY = sqrtOmega0%*%Y
  whiteXaug = scale(sqrtOmega0%*%Xaug)
  ridge1 = cv.glmnet(x = whiteXaug, y = whiteY, alpha = 0)
  ridgeCoef = coef(ridge1, s = ridge1$lambda.min)[-1]

  alasso1 = cv.glmnet(whiteXaug, whiteY, alpha = 1, penalty.factor = 1/abs(ridgeCoef))
  optimalBeta = coef(alasso1, s = "lambda.min")[-1]
  res = list(fit = alasso1, optimalBeta = optimalBeta)
  return(res)
}

#' Run the selection method through SCAD penalty (Fan, 2001)
#'
#' @param Xaug Design matrix of original data column augmented by knockoffs
#' @param Y vectorized response vector
#' @param n the number of samples
#' @param Omega0 the estimated omega matrix
#' @return The result from the SCAD fit
SCADSelection = function(Xaug, Y, n, Omega0 = NULL){
  for(i in 1:n){
    if(i == 1){
      blockOmega0 = Omega0
    }
    else{
      blockOmega0 = bdiag(blockOmega0, Omega0)
    }
  }
  sqrtOmega0 = sqrtm(blockOmega0)

  whiteY = sqrtOmega0%*%Y
  whiteXaug = sqrtOmega0%*%Xaug

  # Run SCAD
  scad.l2 = cv.ncvreg(whiteXaug, whiteY, penalty = "SCAD")
  optimalBeta = coef(scad.l2, s = "lambda.min")[-1]
  scadFit = scad.l2
  res = list(optimalBeta = optimalBeta, scadFit = scadFit)
  return(res)
}
