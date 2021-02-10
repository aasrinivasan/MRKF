#' Run the multiple response knockoff filter
#'
#' @param W OTU count matrix
#' @param Y matrix of responses
#' @param q target FDP
#' @param penalty type of penalty for selection

MRKF = function(W, Ymat, q, penalty = "lasso"){
  ### Load Packages ###
  print("Loading Packages")
  loadPackages()
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

    predictorNames = colnames(Z)[-pstar]

    rownames(selectedMatrix) = predictorNames
    rownames(selectedMatrixp) = predictorNames
    colnames(selectedMatrix) = responseNames
    colnames(selectedMatrixp) = responseNames
    print(selectedMatrix)
    print("Returing Selection Set")
    results = list("S" = selectedMatrix, "Sp"= selectedMatrixp)
    return(results)
  }
}
