module ExternalSolver

  using OSSDP,OSSDPTypes #SCS, JuMP
  export solveWithMOSEK, solveWithSCS

  function solveWithMOSEK(MOSEKobj,MOSEKtime,m,n)

    return  OSSDPTypes.OSSDPResult(zeros(n),zeros(m),zeros(m),zeros(m),MOSEKobj,1,:Solved,MOSEKtime,0.,0.,0.,0.,0.,0.);

  end



# function solvewithSCS(P,q,A,b,Ks,maxIter)
#   m,n = size(A)
#   cDim = Int(sqrt(m))
#   model = Model(solver=SCSSolver(verbose=0,max_iters = maxIter))
#   @variable(model, x[1:n])
#   @variable(model, S[1:cDim,1:cDim],SDP)
#   s = vec(S)
#   @objective(model, Min,q'*x)
#   @constraint(model, A*x.+s .== b)

#   tic()
#   status = JuMP.solve(model)
#   solveTime = toq()

#   cost = getobjectivevalue(model)
#   status == :Optimal ? (status=:Solved) : (status=:Max_iter_reached)

#   return QOCS.Result(getvalue(x),getvalue(s),[0.],[0.],cost,1,status,solveTime,0.,0.,0.,0.,0.,0.);
# end


end