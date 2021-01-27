#' Load all packages needed
#' @description loads all required external packages
#' @keywords internal
#' @noRd
#' @export
loadPackages = function(){
  library(mvtnorm)
  library(knockoff)
  library(QUIC)
  library(R.utils)
  library(MethylCapSig)
  library(GUniFrac)
  library(energy)
  library(compositions)
  library(doParallel)
  library(foreach)
  library(parallel)
  library(ncvreg)
  library(glmnet)
  library(CVglasso)
  library(ncvreg)
  library(Matrix)
  library(matrixcalc)
  library(expm)
  library(parcor)
}
