loadPackages()
set.seed(123)
# Initialize Parameters
n = 100
p = 200
pstar = p + 1
K = 2
A = 3
k = 5
q = .2
sp = .1

# Generate covariance of W
sigma = matrix(0, nrow = pstar, ncol = pstar)
for(i in 1:pstar){
  for(j in 1:pstar){
    sigma[i,j] = .7^(abs(i-j))
  }
}

# generate W
means = rep(1, pstar)
W = mvlognormal(n, means, Sigma = diag(sigma), R = sigma)
colnames(W) = paste("Taxa", 1:ncol(W))

# generate B
means = rep(1, pstar)

print(A)
Aset = seq(-A,A)
if(length(which(Aset==0))!=0){
  Aset = Aset[-which(Aset==0)]
}
B = matrix(0, nrow = p, ncol = K)
totalSig = p*K
nSig = floor(totalSig*sp)
signalEntries = sort(sample(1:totalSig, size = nSig, replace = F))
Bset = sample(Aset, nSig, replace = T)
B[signalEntries] = Bset


# generate error matrix
errorCovariance = matrix(0,K, K)
for(i in 1:K){
  for(j in 1:K){
    errorCovariance[i,j] = errorRho^abs(i-j)
  }
}

E = matrix(0, nrow = n, ncol = K)
for(i in 1:n){
  ei = rnorm(K, mean = rep(0, K), errorCovariance)
  E[i,] = ei
}

# generate Z
Z = acomp(W)
Z = as.matrix(Z)

# generate X
ref = Z[,pstar]
X = matrix(0, nrow = n, ncol = p)
for(t in 1:ncol(X)){
  X[,t] = log(Z[,t]/ref)
}

Y = X%*%B + E
colnames(Y) = paste("Response", 1:ncol(Y))
result = MRKF(W, Y, q, "lasso")
print(result)
