#------------------------------------------------------------------------------
function excise(x)
#notice: if !any(vv), then this does NOT create a copy, so changing the output will
#will change the input

  vv = vec(any(isnan.(x),2))

  if any(vv)              #only keep rows with no NaNs
    vvb = broadcast(!,vv)               #use .~vv in 0.6?
    x = x[vvb,:]                        #use view(x,vvb,:) in the future?
  end

  return x

end
#------------------------------------------------------------------------------
