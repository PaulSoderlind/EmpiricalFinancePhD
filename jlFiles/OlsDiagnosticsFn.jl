function OlsDiagnosticsFn(y,x,u,m=1)
#OlsDiagnosticsFn  Diagnostic tests of OLS residuals
#
#
#
#  Usage:    [AutoCorr,DW,BoxPierce,White,Regr] = OlsDiagnosticsFn(y,x,u,m)
#
#  Input:    y         Tx1, dependent variable
#            x         Txk, regressors
#            u         Tx1, residuals
#            m         scalar, number of lags in autocorrelation and Box-Pierce test
#
#
#  Output:   AutoCorr  mx2, sqrt(T)*autorrelation and p-value
#            DW        1x2, DW statistic and NaN
#            BoxPierce 1x3, Box-Pierce statistic, p-value, df
#            White     1x3, White's test static, p-value, df
#            Regr      1x3, chi2-stat, p-value and df for test of all slope coefficients
#
#
#
#
#  Paul.Soderlind@unisg.ch
#------------------------------------------------------------------------------

  ux = excise([y u x])
  y  = ux[:,1]
  u  = ux[:,2]
  x  = ux[:,3:end]

  T = size(x,1)
  k = size(x,2)

  Stdu = std(u)

  rho = autocor(u,1:m)

  AutoCorrStat = sqrt(T)*rho
  pval         = 2*(1 - cdf(Normal(0,1),abs.(AutoCorrStat)))
  AutoCorr     = [AutoCorrStat pval]

  BPStat       = T*sum(rho.^2,1)
  pval         = 1 - cdf(Chisq(m),BPStat)
  BoxPierce    = [BPStat pval m]

  udiff  = u - [NaN;u[1:end-1]]
  dwStat = mean(excise(udiff).^2,1)/Stdu^2
  pval   = NaN
  DW     = [dwStat pval]

  psi = Array{Float64}(T,0)                #matrix of cross products of x
  for i = 1:k
    for j = i:k
      psi = [psi (x[:,i].*x[:,j])]         #All cross products, incl own
    end
  end

  (b,res,yhat,Covb,R2a,) = OlsFn(u.^2,psi)   #White's test for heteroskedasticity
  WhiteStat = T*R2a./(1-R2a)
  pval      = 1 - cdf(Chisq(size(psi,2)-1),WhiteStat)
  White     = [WhiteStat pval (size(psi,2)-1)]

  (b,res,yhat,Covb,R2a,) = OlsFn(y,x)           #test of regression
  RegrStat = T*R2a./(1-R2a)
  pval     = 1 - cdf(Chisq(size(x,2)-1),RegrStat)
  Regr     = [RegrStat pval (size(x,2)-1)]

  return AutoCorr,DW,BoxPierce,White,Regr

end
#------------------------------------------------------------------------------
