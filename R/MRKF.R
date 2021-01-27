#' Run the multiple response knockoff filter
#'
#' @param W OTU count matrix
#' @param Y matrix of responses
#' @param q target FDP
#' @param penalty type of penalty for selection

multilayerKF = function(W, Y, q, penalty = "lasso"){
  ### Load Packages ###
  loadPackages()

  K = ncol(Y)
  n = nrow(W)
  p = ncol(W)
  pstar = p + 1
  # Compute Z
  Z = acomp(W)
  ref = Z[,pstar]

  # log transformation
  X = matrix(0, nrow = n, ncol = p)
  for(t in 1:ncol(X)){
    X[,t] = log(Z[,t]/ref)
  }

  # Generate Knockoffs
  Xkf = create.second_order(X)
  Im = diag(K)
  for(i in 1:n){
    Xi = kronecker(t(X[i,]), Im)
    Xki = kronecker(t(Xkf[i,]),Im)
    tXi = t(Xi)
    tXki = t(Xki)
    if(i == 1){
      Xmat = tXi
      Xkmat = tXki
    }
    else{
      Xmat = cbind(Xmat, tXi)
      Xkmat = cbind(Xkmat, tXki)
    }
  }
  Xmat = t(Xmat)
  Xkmat = t(Xkmat)

  # Generate augmented
  Xaug = cbind(Xmat, Xkmat)

  # Stage 1
  step1Omega = Step1(Xaug, Y, n, K, Lasso = Lasso)

  # Stage 2
  step2Result = Step2(Xaug, Y, step1Omega, n, Lasso = Lasso)

  if(penalty == "scad"){
    selectStep2 = step2Result$optimalBeta
    Wstat = computeLCD(selectStep2)
  }
  if(penalty == "lasso"){
    selectStep2 = step2Result$optimalBeta
    Wstat = computeLCD(selectStep2)
  }
  if(penalty == "alasso"){
    selectStep2 = step2Result$optimalBeta
    Wstat = computeLCD(selectStep2)
  }

  t = knockoff.threshold(Wstat, fdr = q, offset = 0)
  tp = knockoff.threshold(Wstat, fdr = q, offset = 1)
  S = which(Wstat >= t)
  Sp = which(Wstat >= tp)

  # Return selection sets
  results = list("S" = S, "Sp"=Sp, "t" = t, "tp" = tp, "W" = Wstat)
  return(results)
}
