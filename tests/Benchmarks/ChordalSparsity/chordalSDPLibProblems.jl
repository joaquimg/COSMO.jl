# Test file to test chordal sparsity exploitation of solver on 4 SDPLib problems

workspace()
include("../../../src/Solver.jl")
using OSSDP, Base.Test,Helper


# test solver on the following problems: #maxG32, maxG51, thetaG51, qpG51





K = OSSDPTypes.Cone(Kf,Kl,Kq,Ks)




# solve problem with sparsity turned on
settings1 = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=10000,verbose=false,adaptive_rho=true)
res_sparsity,nothing = OSSDP.solve(P,q,A,b,K,settings1);


# solve problem without sparsity exploitation
settings2 = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=10000,verbose=false,adaptive_rho=true)
res_standard,nothing = OSSDP.solve(P,q,A,b,K,settings2);


# compare results
@test abs(res_sparsity-res_standard) < 1e-3