"""
    OlsnFn(y,x)

LS of y on x; for one dependent variable or SURE with same regressors

# Usage
(b,res,yhat,Covb,R2a,T) = OlsnFn(y,x)

# Input
- `y::Array`:     Txn matrix of the dependent variables
- `x::Array`:     Txk matrix of regressors (including deterministic ones)

# Output
- `b::Array`:     kxn matrix, regression coefficients
- `res::Array`:   Txn matrix, residuals y - yhat
- `yhat::Array`:  Txn matrix, fitted values x*b
- `Covb::Array`:  matrix, covariance matrix of vec(b) =[beq1;beq2;...]
- `R2a::Array`:   1xn vector, R2 values

"""
function OlsnFn(y,x)

  (T,n) = size(y,1,2)

  b     = x\y
  yhat  = x*b
  res   = y - yhat

  Covres = cov(res)*(T-1)/T
  Covb   = kron(Covres,inv(x'x))
  R2a    = 1 .- var(res,1)./var(y,1)

  return b,res,yhat,Covb,R2a

end
