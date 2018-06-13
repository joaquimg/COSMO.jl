workspace()
include("../../../src/Solver.jl")
using OSSDP, Base.Test


function projectCone(x)
  n = Int(sqrt(length(x)))

  # handle 1D case
  if size(x,1) == 1
    return max.(x,0)
  else
    # recreate original matrix from input vectors
    X = reshape(x,n,n)
    X = X./2
    X = X+X'

    # compute eigenvalue decomposition
    F = eigfact(X)

    ind = find(x-> x>0, F[:values])
    Λ = diagm(F[:values])
    UsE = F[:vectors][:,ind]*sqrt.(Λ[ind,ind])
    Xp = UsE*UsE'
    return vec(Xp)
  end
end

rng = MersenneTwister(123123)
maxIter = 10
A = [1. 1 1 1; 1 1 1 0; 1 1 1 1;1 0 1 1 ]
A = vec(A)
m, = size(A)
n = 1

# example doesnt use maximal cliques, therefore aigler's theorem wont work!
Strue = zeros(4,4)
Strue[1:3,1:3] = [6 3 5; 3 6 7; 5 7 14.]
Strue[[1,3,4],[1,3,4]] = [6. 5 1;5 14 1;1 1 3]

xtrue = [1.]
b = A*xtrue[1]+vec(Strue)

P = zeros(1,1)
Ytrue = Helper.generatePosDefMatrix(4,rng)
ytrue = vec(Ytrue)
q = (-P*xtrue -  A'*ytrue)[:]

σ = 1e-6
ρ = 0.1

E1 = [1 0 0 0;0 1 0 0;0 0 1 0] #C1:{1 2 3}
E2 = [1 0 0 0;0 0 1 0;0 0 0 1] #C2: {1 3 4}
H1 = kron(E1,E1)
H2 = kron(E2,E2)
D = H1'*H1 + H2'*H2

# xt = zeros(n,1)
# ν = zeros(m,1)
# x = zeros(n,1)
# s = zeros(m,1)
# st = zeros(m,1)
# μ = zeros(m,1)

# find solution with QOCS
K = OSSDPTypes.Cone(0,0,[],[16])
settings1 = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=10000,verbose=false,adaptive_rho=true)
res,nothing = OSSDP.solve(P,q,A,b,K,settings1);



# ---------------------------------------------
# Only split projection
# ---------------------------------------------

# # remove zeros from A, D to keep system full rank
usedVar = find(x->A[x] != 0., collect(1:m))
Ar = A[usedVar]
br = b[usedVar]
d = diag(D)
Dr = diagm(d[usedVar])
H1r = zeros(9,14)
H1r[1,1] = H1r[2,2] = H1r[3,3] = H1r[4,5] = H1r[5,6] = H1r[6,7] = H1r[7,8] = H1r[8,9] = H1r[9,10] = 1
H2r = zeros(9,14)
H2r[1,1] = H2r[2,3] = H2r[3,4] = H2r[4,8] = H2r[5,10] = H2r[6,11] = H2r[7,12] = H2r[8,13] = H2r[9,14] = 1

numCLiqueNodes = 18
strue = vec(Strue)[usedVar]
nz = length(usedVar)
xt = zeros(n,1)
ν = zeros(nz,1)
x = zeros(n,1)
s = zeros(numCLiqueNodes,1)
st = zeros(numCLiqueNodes,1)
μ = zeros(numCLiqueNodes,1)
s_sum = zeros(nz,1)
μ_sum = zeros(nz,1)
IND = [[1;2;3;5;6;7;8;9;10],[1;3;4;8;10;11;12;13;14] ]
CIND = [IND[1];IND[2]]

# formulate KKT
M = [P+σ*eye(n) Ar'; Ar -1/ρ*Dr]
F = ldltfact(sparse(M))
s1 = 1
s2 = 1
for iii = 1:1500

  RHS = [-q+σ*x;br-s_sum + 1/ρ.*μ_sum]
  sol = F\RHS

  xt = sol[1:n]
  ν = sol[n+1:end]

  x = xt
  st = s - (ν[CIND] + μ)./ρ

  proj = st + μ./ρ
  proj1 = proj[1:9]
  proj2 = proj[10:18]

  s[1:9] = projectCone(proj[1:9])
  s[10:18] = projectCone(proj[10:18])

  s_sum =  zeros(14,1)
  s_sum[IND[1]] = s[1:9]
  s_sum[IND[2]] += s[10:18]

  μ = μ + ρ.*(st - s)

  μ_sum =  zeros(14,1)
  μ_sum[IND[1]] = μ[1:9]
  μ_sum[IND[2]] += μ[10:18]

end

sol1_s = zeros(m)
sol1_s[usedVar] = s_sum
sol1_x = x
cost1 = q'*x
sol1_ν = ν
sol1_μ = μ_sum


# ---------------------------------------------
# No decomposition
# ---------------------------------------------
xt = zeros(n,1)
ν = zeros(m,1)
x = zeros(n,1)
s = zeros(m,1)
st = zeros(m,1)
μ = zeros(m,1)

# formulate KKT
M = [P+σ*eye(n) A'; A -1/ρ*eye(m)]
F = ldltfact(sparse(M))

for iii = 1:1500
  RHS = [-q+σ*x;b-s + 1/ρ.*μ]

  sol = F\RHS

  xt = sol[1:n]
  ν = sol[n+1:end]

  x = xt
  st = s - (ν+μ)./ρ

  proj = st + μ./ρ


  s = projectCone(proj)

  μ = μ + ρ.*(st - s)
end

sol2_s = s
sol2_x = x
cost2 = q'*x
sol2_ν = ν
sol2_μ = μ

# ---------------------------------------------
# Decompose all steps
# ---------------------------------------------
nz = length(usedVar)
xt = zeros(n,1)
ν = zeros(nz,1)
x = zeros(n,1)

s = zeros(nz,1)
s1 = zeros(9,1)
s2 = zeros(9,1)


μ = zeros(nz,1)
μ1 = zeros(9,1)
μ2 = zeros(9,1)


# formulate KKT
M = [P+σ*eye(n) Ar'; Ar -1/ρ*Dr]
F = ldltfact(sparse(M))

for iii = 1:1500

  RHS = [-q+σ*x;br-s + 1/ρ.*μ]
  sol = F\RHS

  xt = sol[1:n]
  ν = sol[n+1:end]

  x = xt
  st1 = s1 - (H1r*ν+μ1)./ρ
  st2 = s2 - (H2r*ν+μ2)./ρ

  proj1 = st1 + μ1./ρ
  proj2 = st2 + μ2./ρ

  s1 = projectCone(proj1)
  s2 = projectCone(proj2)

  s =  H1r'*s1 + H2r'*s2

  μ1 = μ1 + ρ.*(st1 - s1)
  μ2 = μ2 + ρ.*(st2 - s2)


  μ = H1r'*μ1 + H2r'*μ2
end

sol3_s = zeros(m)

sol3_s[usedVar] = s
sol3_x = x
cost3 = q'*x
sol3_ν = ν
sol3_μ = μ
@testset "Simple ADMM loop" begin
  @test abs(res.cost - cost1) < 1e-2
  @test abs(res.cost - cost2) < 1e-2
  @test abs(res.cost - cost3) < 1e-2
  @test minimum(eig(reshape(sol1_s,4,4))[1]) > -1e-5
  @test norm(sol1_s - sol2_s,Inf) < 1e-3
  @test norm(sol1_s - sol3_s,Inf) < 1e-3
  @test norm(sol2_s - sol3_s,Inf) < 1e-3
end