# workspace()
# include("../../../src/Solver.jl")
using Helper, OSSDP, Base.Test


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
A = [1. 1 0; 1 1 1; 0 1 1]
A = vec(A)
m, = size(A)
n = 1

Strue = [5. 1. 0; 1. 4. 1; 0. 1 3]
xtrue = [1.]
b = A*xtrue[1]+vec(Strue)

P = zeros(1,1)
Ytrue = generatePosDefMatrix(3,rng)
ytrue = vec(Ytrue)
q = (-P*xtrue -  A'*ytrue)[:]

σ = 1e-6
ρ = 0.1

E1 = [1 0 0;0 1 0]
E2 = [0 1 0;0 0 1]
H1 = kron(E1,E1)
H2 = kron(E2,E2)
D = H1'*H1 + H2'*H2

xt = zeros(n,1)
ν = zeros(m,1)
x = zeros(n,1)
s = zeros(m,1)
st = zeros(m,1)
μ = zeros(m,1)


# find solution with QOCS
K = OSSDPTypes.Cone(0,0,[],[9])
settings1 = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=10000,verbose=false,adaptive_rho=true)
res,nothing = OSSDP.solve(P,q,A,b,K,settings1);

# formulate KKT
M = [P+σ*eye(n) A'; A -1/ρ*eye(m)]
F = ldltfact(sparse(M))

for iii = 1:500
  RHS = [-q+σ*x;b-s + 1/ρ.*μ]

  sol = F\RHS

  xt = sol[1:n]
  ν = sol[n+1:end]

  x = xt
  st = s - (ν+μ)./ρ
  s = projectCone(st + μ./ρ)

  μ = μ + ρ.*(st - s)
end


@test abs(res.cost - q'*x) < 1e-2