"""
    OlsSureFn(Y,X,m=0)

LS of Y on X; for one n dependent variables, same regressors

# Usage
(b,res,Yhat,Covb,R2) = OlsSureFn(Y,X,m)

# Input
- `Y::Array`:     Txn, the n dependent variables
- `X::Array`:     Txk matrix of regressors (including deterministic ones)
- `m::Int`:       scalar, bandwidth in Newey-West

# Output
- `b::Array`:     n*kx1, regression coefficients
- `u::Array`:     Txn, residuals Y - Yhat
- `Yhat::Array`:  Txn, fitted values X*b
- `V::Array`:     matrix, covariance matrix of vec(b)
- `R2::Number`:  n vector, R2 value

"""
function OlsSureFn(Y,X,m=0)

    (T,n) = (size(Y,1),size(Y,2))
    k     = size(X,2)

    b     = X\Y
    Yhat  = X*b
    u     = Y - Yhat

    g     = zeros(T,n*k)
    for i = 1:n
      vv      = (1+(i-1)*k):(i*k)   #1:k,(1+k):2k,...
      g[:,vv] = X.*u[:,i]           #moment conditions for Y[:,i] regression
    end

    S0    = NWFn(g,m)            #Newey-West covariance matrix
    Sxx   = X'X/T
    Sxx_1 = kron(I(n),inv(Sxx))
    V     = Sxx_1 * S0 * Sxx_1/T

    R2   = 1 .- var(u,dims=1)./var(Y,dims=1)

    return b, u, Yhat, V, R2

end
