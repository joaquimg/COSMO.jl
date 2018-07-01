# Unit test to check if the function findStackingMatrix() works properly
workspace()
include("../../../src/Types.jl")
include("../../../src/Graph.jl")
include("../../../src/Tree.jl")
include("../../../src/ChordalSparsity.jl")


using Base.Test,ChordalSparsity







rng = MersenneTwister(7821)

Kf = 1
Kl = 2
Kq = [3]
Ks = [9 16]
K = OSSDPTypes.Cone(Kf,Kl,Kq,Ks)

cliqueSets = Array{ChordalSparsity.CliqueSet,1}(2)
cliqueSets[1] = CliqueSet([[1;2],[2;3]])
cliqueSets[2] = CliqueSet([[1;2;3],[2;3;4]])

 H,Ks2 = ChordalSparsity.findStackingMatrix(K,cliqueSets)

sbar =[ [1;1;2;1;2;3];[1;2;4;5;5;6;8;9];[1;2;3;5;6;7;9;10;11;6;7;8;10;11;12;14;15;16]]

s = H*sbar