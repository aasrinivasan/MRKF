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

# W covariance
sigma = matrix(0, nrow = pstar, ncol = pstar)

for(i in 1:pstar){
  for(j in 1:pstar){
    sigma[i,j] = .5^(abs(i-j))
  }
}

means = rep(1, pstar)
W = mvlognormal(n, means, Sigma = diag(sigma), R = sigma)
colnames(W) = paste("Taxa", 1:ncol(W))
# generate B
B = rep(0,p)
signals = sample(1:p,k)
signalValue = runif(k,min = -A, max = A)
B[signals] = signalValue

# gererate error
eps = rnorm(n, 0, 1)
eps2 = rnorm(n,0,1)

# generate Z
Z = acomp(W)
Z = as.matrix(Z)

# generate X
ref = Z[,pstar]
X = matrix(0, nrow = n, ncol = p)
for(t in 1:ncol(X)){
  X[,t] = log(Z[,t]/ref)
}

# generate Y
Y = X%*%B + eps
Y2 =  X%*%B + eps

Ymatrix = cbind(Y, Y2)
colnames(Ymatrix) = paste("Response", 1:ncol(Ymatrix))
result = MRKF(W, Ymatrix, q, "alasso")
print(result)
