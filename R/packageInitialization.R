#' Load all packages needed
#' @description loads all required external packages
#' @keywords internal
#' @noRd
#' @export packageInitialization
packageInitialization = function(){
  library(mvtnorm)
  library(knockoff)
  library(MethylCapSig)
  library(GUniFrac)
  library(energy)
  library(compositions)
  library(doParallel)
  library(foreach)
  library(ncvreg)
  library(glmnet)
  library(CVglasso)
  library(ncvreg)
  library(Matrix)
  library(matrixcalc)
  library(expm)
}
