#------------------------------------------------------------------------------
"""
Ols2Fn(y,x,m=0)

LS of y on x for one dependent variable or SURE with same regressors, gives a NW
covariance matrix.

 # Input
- `y::Array`:  Tx1 or Txn matrix of the dependent variables
- `x::Array`:  Txk matrix of regressors (including deterministic ones)
- `m::Int`:    scalar, bandwidth in the Newey-West covariance matrix
               Use m = 0 to get White's covariance matrix

 # Output
- `b::Array`:       kxn matrix, regression coefficients
- `u::Array`:       Tx1 or Txn matrix, residuals y - yhat
- `yhat::Array`:    Tx1 or Txn matrix, fitted values x*b
- `CovbLS::Array`:  matrix, covariance matrix of vec(b) = [beq1;beq2...]
- `R2a::Array`:     1xn vector, R2 values
- `T::Int`:         scalar, number of obs
- `CovbNW::Array`:  matrix, Newey-West covariance matrix. Use m=0 for White's

# Requires
- HDirProdFn
- NWFn

"""
#----------------------------------------------------------------------------

function Ols2Fn(y,x,m=0)

  T = size(y,1)
  n = size(y,2)

  b    = x\y
  yhat = x*b
  u    = y - yhat

  Covres = cov(u)*(T-1)/T
  CovbLS = kron(Covres,inv(x'x))
  R2a    = 1 - var(u,1)./var(y,1)

  #g      = x .* u                  #x.*u
  g      = HDirProdFn(u,x)
  S      = NWFn(g,m)              #Newey-West
  D      = -x'x/T
  CovbNW = inv(D'inv(S)*D)/T     #Cov(b)

  return b,u,yhat,CovbLS,R2a,T,CovbNW

end
#------------------------------------------------------------------------------
