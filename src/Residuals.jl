module Residuals
using OSSDPTypes
export calculateResiduals, maxResComponentNorm, hasConverged,investigateMaxNorms!, calculateOriginalResiduals!

function calculateResiduals(ws::OSSDPTypes.WorkSpace,settings::OSSDPTypes.OSSDPSettings, IGNORESCALING_FLAG::Bool=false)
  if (settings.scaling != 0 && !IGNORESCALING_FLAG)
    r_prim = norm(ws.sm.Einv*(ws.p.A*ws.x+ws.s-ws.p.b),Inf)
    # ∇f0 + ∑ νi ∇hi(x*) == 0 condition
    r_dual = norm(ws.sm.cinv*ws.sm.Dinv*(ws.p.P*ws.x + ws.p.q -ws.p.A'*ws.μ),Inf)
  end
  if (settings.scaling == 0 || IGNORESCALING_FLAG )
    r_prim = norm(ws.p.A*ws.x+ws.s-ws.p.b,Inf)
    r_dual = norm(ws.p.P*ws.x + ws.p.q -ws.p.A'*ws.μ,Inf)
  end
  # FIXME: Why is it -A'μ ?
  return r_prim,r_dual
end




function calculateOriginalResiduals!(ws::OSSDPTypes.WorkSpace,residualInfo)
  mO = ws.ci.originalM
  nO = ws.ci.originalN

  A = (ws.sm.Einv*ws.p.A*ws.sm.Dinv)[1:mO,1:nO]
  P = ws.sm.cinv*(ws.sm.Dinv*ws.p.P*ws.sm.Dinv)[1:nO,1:nO]
  q = ws.sm.cinv*(ws.sm.Dinv*ws.p.q)[1:nO]
  b = (ws.sm.Einv*ws.p.b)[1:mO]

  # calculate residuals
  residualInfo.rPrim = norm(A*ws.x + ws.s - b,Inf)
  residualInfo.rDual = norm(P*ws.x  + q - A'*ws.μ,Inf)


  # calculate max terms
  maxNormPrim1 = norm(A*ws.x,Inf)
  maxNormPrim2 = norm(ws.s,Inf)
  maxNormPrim3 = norm(b,Inf)
  maxNormDual1 = norm(P*ws.x,Inf)
  maxNormDual2 = norm(q,Inf)
  maxNormDual3 = norm(A'*ws.μ,Inf)

  residualInfo.rPrimMax = [maxNormPrim1;maxNormPrim2;maxNormPrim3]
  residualInfo.rDualMax = [maxNormDual1;maxNormDual2;maxNormDual3]

    return nothing
end


function maxResComponentNorm(ws::OSSDPTypes.WorkSpace,settings::OSSDPTypes.OSSDPSettings, IGNORESCALING_FLAG::Bool=false)
  if (settings.scaling != 0 && !IGNORESCALING_FLAG)
    maxNormPrim = max.(norm(ws.sm.Einv*ws.p.A*ws.x,Inf),  norm(ws.sm.Einv*ws.s,Inf), norm(ws.sm.Einv*ws.p.b,Inf))
    maxNormDual = max.(norm(ws.sm.cinv*ws.sm.Dinv*ws.p.P*ws.x,Inf), norm(ws.sm.cinv*ws.sm.Dinv*ws.p.q,Inf), norm(ws.sm.cinv*ws.sm.Dinv*ws.p.A'*ws.μ,Inf) )
  end
  if (settings.scaling == 0 || IGNORESCALING_FLAG)
    maxNormPrim = max.(norm(ws.p.A*ws.x,Inf),norm(ws.s,Inf), norm(ws.p.b,Inf))
    maxNormDual = max.(norm(ws.p.P*ws.x,Inf), norm(ws.p.q,Inf),norm(ws.p.A'*ws.μ,Inf) )
  end
  return maxNormPrim, maxNormDual
end

function investigateMaxNorms!(ws,settings,residualInfo)
  IGNORESCALING_FLAG = false
  maxNormPrim1 = maxNormPrim2 = maxNormPrim3 = maxNormDual1 = maxNormDual2 = maxNormDual3 = 0.
  if (settings.scaling != 0 && !IGNORESCALING_FLAG)
    maxNormPrim1 = norm(ws.sm.Einv*ws.p.A*ws.x,Inf)
    maxNormPrim2 = norm(ws.sm.Einv*ws.s,Inf)
    maxNormPrim3 = norm(ws.sm.Einv*ws.p.b,Inf)
    maxNormDual1 = norm(ws.sm.cinv*ws.sm.Dinv*ws.p.P*ws.x,Inf)
    maxNormDual2 = norm(ws.sm.cinv*ws.sm.Dinv*ws.p.q,Inf)
    maxNormDual3 = norm(ws.sm.cinv*ws.sm.Dinv*ws.p.A'*ws.μ,Inf)
  end
  if (settings.scaling == 0 || IGNORESCALING_FLAG)
    maxNormPrim1 = norm(ws.p.A*ws.x,Inf)
    maxNormPrim2 = norm(ws.s,Inf)
    maxNormPrim3 = norm(ws.p.b,Inf)
    maxNormDual1 = norm(ws.p.P*ws.x,Inf)
    maxNormDual2 = norm(ws.p.q,Inf)
    maxNormDual3 = norm(ws.p.A'*ws.μ,Inf)
  end
  residualInfo.rPrimMax_aug = [maxNormPrim1;maxNormPrim2;maxNormPrim3]
  residualInfo.rDualMax_aug = [maxNormDual1;maxNormDual2;maxNormDual3]
  return nothing
end



function hasConverged(ws::OSSDPTypes.WorkSpace,settings::OSSDPTypes.OSSDPSettings,r_prim::Float64,r_dual::Float64)
  maxNormPrim, maxNormDual = maxResComponentNorm(ws,settings)
  ϵ_prim = settings.eps_abs + settings.eps_rel * maxNormPrim
  ϵ_dual = settings.eps_abs + settings.eps_rel * maxNormDual

  # if an optimal objective value was specified for the problem check if current solution is within specified accuracy
  objTrueFLAG = true
  if !isnan(settings.objTrue)
    currentCost = ws.sm.cinv*(1/2 * ws.x'*ws.p.P*ws.x + ws.p.q'*ws.x)[1]
    objTrueFLAG = abs(settings.objTrue-currentCost) <= settings.objTrueTOL
  end

  return ( r_prim < ϵ_prim  && r_dual < ϵ_dual && objTrueFLAG)
end

end #MODULE