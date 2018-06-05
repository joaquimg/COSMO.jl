
workspace()
include("../../../src/Solver.jl")
include("SparseSDPs.jl")
using OSSDP,SparseSDPs
using Base.Test
rng = MersenneTwister(3232)


# generate a random test problem with banded sdp sparsity pattern
m = 3
K = OSSDPTypes.Cone(2,3,[3 5],[8 10])
bandwidth = [2 3]
density = 0.45
P,q,A,b,K = generateBandedSDP(rng,m,K,bandwidth,-20.,20.,density);






# # solve problem with sparsity turned on
settings1 = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=10000,verbose=false,adaptive_rho=true)
res_sparsity,nothing = OSSDP.solve(P,q,A,b,K,settings1);
nothing

# # solve problem without sparsity exploitation
# settings2 = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=10000,verbose=false,adaptive_rho=true)
# res_standard,nothing = OSSDP.solve(P,q,A,b,K,settings2);


# # compare results
# @test abs(res_sparsity-res_standard) < 1e-3


