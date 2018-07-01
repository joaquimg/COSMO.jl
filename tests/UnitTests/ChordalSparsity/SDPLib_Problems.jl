workspace()
include("../../../src/Solver.jl")
using JLD
using Base.Test
using OSSDP
using ChordalSparsity, OSSDPTypes

data = JLD.load("maxG32_inSolverFormat.jld")
A = data["A"]
b = data["b"]
q = data["q"]
P = data["P"]
optVal = data["optVal"]

println("Problem data loaded!")
m,n = size(A)

Kf = 0
Kl = 0
Kq = []
Ks = [m]

# define cone membership
K = Cone(Kf,Kl,Kq,Ks)
settings = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=3,verbose=false,adaptive_rho=true,scaling=10,decompose=true)



A = sparse(rand(16,4))
b = rand(16,1)
A[4,:] = 0
A[14,:] = 0
A[8,:] = 0
b[14] = 0
b[8] = 0
A = dropzeros(A)
tic()
csp = ChordalSparsity.findCommonSparsity(A,b)
sp = GraphModule.Graph(csp,Int(sqrt(size(A,1))),false)
toq()

# chordalInfo = OSSDPTypes.ChordalInfo(size(A,1),size(A,2),K)
# println("Start Chordal Decomposition!")
# Pa, qa, Aa, ba, Ka = ChordalSparsity.chordalDecomposition!(P,q,A,b,K,settings,chordalInfo)
# println("Problem Number of cones increased from $(length(K.s)) to $(length(Ka.s)).")

