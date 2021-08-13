#' Multiple Response Knockoff Filter
#'
#' Run the MRCKF procedure
#'
#' @param W OTU count matrix of dimension (n x p^*)
#' @param Y matrix of responses (n x K)
#' @param q target FDP (0,1)
#' @param penalty type of penalty for selection: "lasso", "alasso", or "scad"
#' @return A list containing the following components:
#' \item{S}{The estimated selection matrix for B under the knockoff threshold}
#' \item{Sp}{The estimated selection matrix for B under the knockoff+ threshold}
#'
#' @details The MRCKF method relies on the model-X formulation which assumes that the underlying
#' alr transformed data \eqn{X} follows a multivariate normal distribution. The model-X formulation generates
#' knockoff copies by estimating the underlying mean and covariance structure. Note that the MRCKF
#' method requires the input of the raw, non-compositional data \eqn{W}.
#'
#' If the number of responses, \eqn{K=1}, then the MRCKF method defaults to a variation of the
#' compostional knockoff filter (CKF) by Srinivasan et al., 2020.
#'
#' The knockoff importance statistic used in this analysis is the lasso coefficient difference statistic.
#'
#' @references
#'
#' Srinivasan, A, Xue, L, Zhan, X. Compositional knockoff filter for high-dimensional regression analysis of microbiome data.
#' Biometrics. 2020; 1– 12. https://doi.org/10.1111/biom.13336
#' @import mvtnorm
#' @import knockoff
#' @import R.utils
#' @import MethylCapSig
#' @import GUniFrac
#' @import energy
#' @import compositions
#' @import ncvreg
#' @import glmnet
#' @import CVglasso
#' @import Matrix
#' @import matrixcalc
#' @import expm
#' @import parcor
#' @import doParallel
#' @import foreach
#' @export MRCKF
MRCKF = function(W, Ymat, q, penalty = "lasso"){
  ### Load Packages ###
  #print("Loading Packages")
  #packageInitialization()
  K = ncol(Ymat)
  n = nrow(W)
  pstar = ncol(W)
  Y = as.vector(Ymat)
  p = pstar - 1

  predictorNames = colnames(W)

  print("Transforming to Compositional Data")
  # Compute Z
  Z = acomp(W)
  Z = as.matrix(Z)
  res = Z[,pstar]

  print("Applying alr Transformation")
  # log transformation
  X = matrix(0, nrow = n, ncol = p)
  for(t in 1:ncol(X)){
    X[,t] = log(Z[,t]/ref)
  }

  if(K == 1){
    print("Running Single Response Setting - Defaulting to Standard Knockoffs")
    result = runSingleKF(X, Y, q, penalty)

    results = list("S" = predictorNames[result$S], "Sp"= predictorNames[result$S])
    return(results)
  }
  else{
    responseNames = colnames(Ymat)
    # Generate Knockoffs
    print("Running Multiple Response Setting")
    print("Generating Knockoffs")
    Xkf = create.second_order(X)
    id = diag(K)
    Xmat = kronecker(id, X)
    Xkmat = kronecker(id, Xkf)
    Xaug = cbind(Xmat, Xkmat)

    # Stage 1
    print("Regression Stage 1")
    step1Omega = Step1(Xmat, Y, n, K, method = penalty)

    # Stage 2
    print("Regression Stage 2")
    sqrtOmega0 = sqrtm(step1Omega)
    blockOmega0 = kronecker(sqrtOmega0, diag(n))
    whiteY = blockOmega0%*%Y
    whiteXaug = blockOmega0%*%Xaug

    step2Result = Step2(whiteXaug, whiteY, step1Omega, n, method = penalty)

    print("Computing Knockoff Statistics")

    selectStep2 = step2Result$optimalBeta
    Wstat = computeLCD(selectStep2)


    t = knockoff.threshold(Wstat, fdr = q, offset = 0)
    tp = knockoff.threshold(Wstat, fdr = q, offset = 1)
    S = which(Wstat >= t)
    Sp = which(Wstat >= tp)

    # Return selection sets
    selected = rep(0, length(Wstat))
    selectedp = rep(0, length(Wstat))

    selected[S] = 1
    selectedp[Sp] = 1

    selectedMatrix = matrix(selected, nrow = p, ncol = K, byrow = F)
    selectedMatrixp = matrix(selectedp, nrow = p, ncol = K, byrow = F)
    print(p)
    print(K)
    predictorNames = colnames(Z)[-pstar]

    rownames(selectedMatrix) = predictorNames
    rownames(selectedMatrixp) = predictorNames
    colnames(selectedMatrix) = responseNames
    colnames(selectedMatrixp) = responseNames
    print("Returning Selection Set")
    results = list("S" = selectedMatrix, "Sp"= selectedMatrixp)
    return(results)
  }
}
