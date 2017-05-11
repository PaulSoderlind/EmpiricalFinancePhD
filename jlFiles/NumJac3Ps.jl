function NumJac3Ps(fun::Function,b0,h=[],Method=99)
#NumJac3Ps    Numerical Jacobian
#
#
#
#
#  Usage:    Jac = NumJac3Ps(fun,b0[,h[,Method]])
#
#  Input:    fun         function to differentiate, generates nx1 output
#            b0          kx1 vector of paramaters, fun(b)
#            h           (optional) scalar or kx1 vector, pertubation of b,
#                        [] gives a good default
#            Method      (optional) 1: forward derivatives; 2: backward derivatives
#                        else: central derivatives
#
#  Output:   Jac         nxk matrix of derivatives, Jac[i,j] = df[i]/db[j]
#
#
#
#  Notice: fun(b0) must be well defined. For instance, it must work even if b0 is
#          an 1-element array
#
#
#  Paul.Soderlind@unisg.ch   2007, to Julia Nov 2015
#------------------------------------------------------------------------------

  eps64 = eps(Float64)

  if isempty(h)
    if any(Method == [1,2])
      h   = sqrt.(eps64)*max.(abs.(b0),1)
      bdown = b0 - h
      bup   = b0 + h
      h     = bup - b0
    else
      h     = eps64^(1/3)*max.(abs.(b0),1)
      bdown = b0 - h
      bup   = b0 + h
      hh    = bup - bdown
    end
  else
    bdown = b0 - h
    bup   = b0 + h
    hh    = bup - bdown
  end


  b0 = vec(collect(b0))                       #-> column vector, for bminus[i]
  k  = length(b0)                             #no. of parameters in fun(b)

  f0 = fun(b0)                               #value of function at b0
  n = length(f0)                             #no. of outputs of fun


  Jac = fill(NaN,(n,k))
  for i = 1:k                               #loop over parameters
    bminus    = deepcopy(b0)
    bminus[i] = bdown[i]
    bplus     = deepcopy(b0)
    bplus[i]  = bup[i]
    if Method == 1
      fplus    = fun(bplus)
      Jac[:,i] = (fplus-f0)/h[i]
    elseif Method == 2
      fminus   = fun(bminus)
      Jac[:,i] = (f0-fminus)/h[i]
    else
      fplus    = fun(bplus)
      fminus   = fun(bminus)
      Jac[:,i] = (fplus-fminus)/hh[i]
    end
  end   #i

  if length(Jac) == 1                   #convert to scalar
    Jac = Jac[1]
  end

  return Jac

end
#------------------------------------------------------------------------------
