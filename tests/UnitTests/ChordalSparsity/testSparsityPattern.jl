
workspace()
include("../../../src/Solver.jl")
include("SparseSDPs.jl")
using OSSDP,SparseSDPs, ChordalSparsity
using Base.Test
rng = MersenneTwister(3232)


# generate a random test problem with banded sdp sparsity pattern
m = 10
K = OSSDPTypes.Cone(0,0,[],[10])
bandwidth = [3]
density = 0.45
P,q,A,b,K = generateBandedSDP(rng,m,K,bandwidth,-20.,20.,density);


As = sparse(full(reshape(A[:,1],10,10)))


As2 = sparse([1 1 0;1 1 1; 0 1 1.])
sp = SparsityPattern(As);

nothing
