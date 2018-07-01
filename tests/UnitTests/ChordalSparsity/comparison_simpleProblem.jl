workspace()
include("../../../src/Solver.jl")
using OSSDP, OSSDPTypes, Base.Test


rng = MersenneTwister(123123)
maxIter = 10
A = [1. 1 1 1; 1 1 1 0; 1 1 1 1;1 0 1 1 ]
A = vec(A)[:,:]
A = sparse(A)
m, = size(A)
n = 1

Strue = zeros(4,4)
Strue[1:3,1:3] = [6 3 5; 3 6 7; 5 7 14.]
Strue[[1,3,4],[1,3,4]] = [6. 5 1;5 14 1;1 1 3]

xtrue = [1.]
b = A*xtrue[1]+vec(Strue)
b = b[:]
P = spzeros(1,1)
Ytrue = Helper.generatePosDefMatrix(4,rng)
ytrue = vec(Ytrue)
q = (-P*xtrue -  A'*ytrue)[:]



# solve with decomposition
K = OSSDPTypes.Cone(0,0,[],[16])
settings = OSSDPSettings(rho=0.1,sigma=1e-6,scaling=0,alpha=1.6,max_iter=2000,verbose=true,adaptive_rho=true,decompose=true)
res_decomposed,nothing = OSSDP.solve(P,q,A,b,K,settings);

# solve once without decomposition
K = OSSDPTypes.Cone(0,0,[],[16])
settings = OSSDPSettings(rho=0.1,sigma=1e-6,scaling=0,alpha=1.6,max_iter=2000,verbose=true,adaptive_rho=true,decompose=false)
res,nothing = OSSDP.solve(P,q,A,b,K,settings);


@testset "Decomposition on simple problem" begin
    @test abs(res.cost-res_decomposed.cost) <= 1e-2
    @test norm(res.x-res_decomposed.x,Inf) <= 1e-2
    @test norm(res.s-res_decomposed.s,Inf) <= 1e-2
    @test minimum(eig(reshape(res.s,4,4))[1]) >= -1e-5
    @test minimum(eig(reshape(res_decomposed.s,4,4))[1]) >= -1e-5
end



