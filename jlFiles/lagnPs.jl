function lagnPs(x,n=1)
#lagnPs   Creates a matrix or vector of lagged values.
#
#
#
#  Usage:    z    = lagnPs(x,n) or
#                 = lagnPs(x)
#
#  Input:    x     Txk matrix
#            n     scalar, order of lag. For instance, 2 puts x(t-2)
#                  on row t; -3 puts x(t+3) on row t.
#
#  Output:   z     Txk matrix of lags
#
#
#
#
#
#  Paul.Soderlind@unisg.ch, to Julia Nov 2015
#--------------------------------------------------------------------------

  (T,k) = size(x,1,2)

  if abs(n) >= T
    error("too many lags or leads")
  end

  missings = fill(NaN,(abs(n),k))
  if n < 0                                    #leads
    z = [ x[1-n:T,:]; missings ]
  elseif n == 0                               #no lag or lead
    z = deepcopy(x)
  elseif n > 0                                #lags
    z = [ missings; x[1:T-n,:] ]
  end

  return z

end
#--------------------------------------------------------------------------
