# Test file to check the kkt solvers
using Test, LinearAlgebra, SparseArrays, Random, COSMO

rng = Random.MersenneTwister(1234)

m = 3
n = 4

function make_test_kkt(P,A,sigma,rho)

    R = length(rho)   == 1 ? ((1.)./rho)*I : Diagonal((1.)./rho)
    S = length(sigma) == 1 ? (sigma)*I : Diagonal(sigma)

    #compute the full KKT matrix
    K    = [P+S A'; A -R]
    return K
end

solver_types = [COSMO.CholmodKKTSolver COSMO.PardisoKKTSolver COSMO.QDLDLKKTSolver]

@testset "$t : KKT solver test" for t in solver_types

    P = sprandn(n,n,0.2)
    P = P'*P
    A = sprandn(m,n,0.5)
    sigma = 1e-5
    rho   = 0.1
    rhs   = randn(m+n)


    for rho in [rand(1)[], rand(m)],
        sigma in [rand(1)[], rand(n)]

        #compute the full KKT matrix
        K = make_test_kkt(P,A,sigma,rho)

        solver = t(P,A,sigma,rho)
        lhs    = copy(rhs)

        #test one time solution
        COSMO.ldiv!(solver,lhs)
        @test (norm(K\rhs - lhs) < 1e-10)

    end

    #test rho update
    rhs = randn(m+n)
    solver = t(P,A,sigma,rho)
    for rho_new in [rand(1)[], rand(m)]

        COSMO.update_rho!(solver,rho_new)
        #compute the full KKT matrix with new rho
        K = make_test_kkt(P,A,sigma,rho_new)
        lhs    = copy(rhs)

        #test updated solution
        ldiv!(solver,lhs)
        @test (norm(K\rhs - lhs) < 1e-10)
    end
end


nothing
