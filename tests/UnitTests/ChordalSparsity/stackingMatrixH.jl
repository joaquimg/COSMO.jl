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

# Tests
 function singleRows(A)
  for iii=1:size(A,1)
    if sum(A[iii,:]) != 1.
      return false
    end
  end
  return true
 end

 function zeroCols(A)
  n = 0
  for jjj=1:size(A,2)
    nnz(A[:,jjj]) == 0 && (n+=1)
  end
  return n
end


 @testset "Stacking Matrix" begin
   @test nnz(H) == size(H,1)
   @test singleRows(H)
   @test zeroCols(H) == (sum(Ks) - sum(Ks2))
end
