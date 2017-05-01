function NWFn(g,m)
 #NWFn    Calculates covariance matrix of sqrt(T)*sample average.

  T = size(g,1)                     #g is Txq
  m = min(m,T-1)                    #number of lags

  g = g .- mean(g,1)                #Normalizing to Eg=0

  S = g'g/T                         #(qxT)*(Txq)
  for s = 1:m
    Gamma_s = g[s+1:T,:]'g[1:T-s,:]/T   #same as Sum[g(t)*g(t-s)',t=s+1,T]
    S       = S  +  ( 1 - s/(m+1) ) * (Gamma_s + Gamma_s')
  end

  return S

end