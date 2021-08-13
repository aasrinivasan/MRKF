### Knockoff Functions
#' @export computeLCD
computeLCD = function(B){
  p = length(B)
  orig = 1:(p/2)
  kf = (p/2+1):p

  oB = abs(B[orig])
  kB = abs(B[kf])
  W = oB-kB
  return(W)
}
#' @export selectFeatures
selectFeatures = function(W, thresh, offset){
  t = knockoff.threshold(W, thresh, offset)
  S = which(W >= t)
  return(S)
}


