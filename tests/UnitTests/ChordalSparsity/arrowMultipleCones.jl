# Test file to compare result for solver when sparsity is exploited vs. standard method

# workspace()
# include("../../../src/Solver.jl")
# include("../../../../../SDP_Benchmark_Problems/Code/SparseSDPs.jl")
# using OSSDP,SparseSDPs
# using Base.Test
# rng = MersenneTwister(9878812)


# generate a random test problem with block-arrow sparsity pattern (similar to CDCS test problem)
m = 3
numCones = 2
nBlk = [2 20]
BlkSize = [2 5]
ArrowWidth = [1 5]
NONZERO_P_FLAG = true
P,q,A,b,Ks = generateArrowMultCones(rng, m,numCones,nBlk,BlkSize,ArrowWidth,NONZERO_P_FLAG);
K = OSSDPTypes.Cone(0,0,[],Ks)






# # solve problem with sparsity turned on
settings1 = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=10000,verbose=true,adaptive_rho=true)
res_sparsity,nothing = OSSDP.solve(P,q,A,b,K,settings1);
nothing

# # solve problem without sparsity exploitation
settings2 = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=10000,verbose=true,adaptive_rho=true)
res_standard,nothing = OSSDP.solve(P,q,A,b,K,settings2);


# # compare results
@test abs(res_sparsity.cost-res_standard.cost) < 1e-3
