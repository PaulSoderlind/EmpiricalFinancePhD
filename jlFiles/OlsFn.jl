#----------------------------------------------------------------------------
"""
    OlsFn(y,x)

LS of y on x; for one dependent variable or SURE with same regressors

# Usage
(b,res,yhat,Covb,R2a,T) = OlsFn(y,x)

# Input
- `y::Array`:     Tx1 or Txn matrix of the dependent variables
- `x::Array`:     Txk matrix of regressors (including deterministic ones)

# Output
- `b::Array`:     kxn matrix, regression coefficients
- `res::Array`:   Tx1 or Txn matrix, residuals y - yhat
- `yhat::Array`:  Tx1 or Txn matrix, fitted values x*b
- `Covb::Array`:  matrix, covariance matrix of vec(b) =[beq1;beq2;...]
- `R2a::Array`:   1xn vector, R2 values
- `T::Int`:       scalar, number of obs

"""
function OlsFn(y,x)

  T = size(y,1)
  n = size(y,2)

  b     = x\y
  yhat  = x*b
  res   = y - yhat

  Covres = cov(res)*(T-1)/T
  Covb   = kron(Covres,inv(x'x))
  R2a    = 1 - var(res,1)./var(y,1)

  return b,res,yhat,Covb,R2a,T

end
#------------------------------------------------------------------------------
