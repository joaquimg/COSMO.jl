# Test file to compare result for solver when sparsity is exploited vs. standard method for a very simple test case

workspace()
include("../../../src/Solver.jl")
using OSSDP, Base.Test,Helper
rng = MersenneTwister(9878812)

# generate a simple 3x3 problem with 2 cones

  A1 = rand(rng,4,4)
  A1[1,3] = A1[1,4] = A1[3,1] = A1[4,1] = 0
  A = vec(A1)
  m, = size(A)
  xtrue = randn(rng)

  S1 = generatePosDefMatrix(2,rng)
  S2 = generatePosDefMatrix(2,rng)
  Strue = blkdiag(sparse(S1),sparse(S2))
  b = A*xtrue+vec(Strue)

  P = sprand(rng,1,1,1.)*(10-0.1)+0.1
  ytrue = rand(rng,m)
  q = full(-P*xtrue - A'*ytrue)[:]
  Kf = 0
  Kl = 0
  Kq = []
  Ks = [m]



K = OSSDPTypes.Cone(Kf,Kl,Kq,Ks)




# solve problem with sparsity turned on
settings1 = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=10000,verbose=false,adaptive_rho=true)
res_sparsity,nothing = OSSDP.solve(P,q,A,b,K,settings1);


# solve problem without sparsity exploitation
settings2 = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=10000,verbose=false,adaptive_rho=true)
res_standard,nothing = OSSDP.solve(P,q,A,b,K,settings2);


# compare results
@test abs(res_sparsity-res_standard) < 1e-3