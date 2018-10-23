using QDLDL, AMD, SparseArrays, Pardiso
import LinearAlgebra: ldiv!

# -------------------------------------
# abstract type defs
# -------------------------------------
abstract type AbstractKKTSolver end

# NB: all concrete examples of this type should
# implement ldiv!, update_rho!, restart!, and should
# implement a construct taking (P,A,sigma,rho) as
# arguments

# -------------------------------------
# some internal utility functions
# -------------------------------------

function _kktutils_check_dims(P,A,sigma,rho)

    n = size(P,1)
    m = size(A,1)

    size(A,2) == n || throw(DimensionMismatch())

    length(rho)   == m || length(rho)   == 1 || throw(DimensionMismatch())
    length(sigma) == n || length(sigma) == 1 || throw(DimensionMismatch())

    return m, n
end

function _kktutils_make_kkt(P,A,sigma,rho,shape::Symbol=:F)

    R = length(rho)   == 1 ? ((1.)./rho)*I : Diagonal((1.)./rho)
    S = length(sigma) == 1 ? (sigma)*I : Diagonal(sigma)

    n = size(P,1)
    m = size(A,1)

    if     shape == :F
        #compute the full KKT matrix
        K = [P+S A'; A -R]

    elseif shape == :U
        #upper triangular
        K = [triu(P)+S  A'; spzeros(eltype(A),m,n)  -R]

    elseif shape == :L
        #lower triangular
        K = [tril(P)+S  spzeros(eltype(A),n,m); A  -R]

    else
        error("Bad matrix shape description")
    end

    return K

end


# -------------------------------------
# QDLDL solver
# -------------------------------------
mutable struct QDLDLKKTSolver <: AbstractKKTSolver

    fact::QDLDL.QDLDLFactorisation
    K::SparseMatrixCSC
    m::Integer
    n::Integer

    function QDLDLKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC,sigma,rho)

        m,n = _kktutils_check_dims(P,A,sigma,rho)
        #NB: qdldl uses triu internally, but it reorders
        #with AMD first.  This way is memory inefficient
        #but saves having to permute the rhs/lhs each
        #time we solve.
        K   = _kktutils_make_kkt(P,A,sigma,rho,:F)
        fact = qdldl(K)

        return new(fact,K,m,n)

    end
end

ldiv!(s::QDLDLKKTSolver, rhs::Vector) = QDLDL.solve!(s.fact,rhs)

restart!(s::QDLDLKKTSolver) = nothing #direct method, nothing to restart

function update_rho!(s::QDLDLKKTSolver, rho::Union{Vector,AbstractFloat})

    didx = diagind(s.K, 0)
    @views s.K[didx[(s.n+1):(s.m+s.n)]] .= -1. ./ rho
    s.fact = qdldl(s.K) #just restart for now.  Very wasteful!
end


# -------------------------------------
# Julia Native solver (CHOLMOD based)
# -------------------------------------
mutable struct CholmodKKTSolver <: AbstractKKTSolver

    fact::SuiteSparse.CHOLMOD.Factor
    K::SparseMatrixCSC
    m::Integer
    n::Integer

    function CholmodKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC,sigma,rho)

        m,n  = _kktutils_check_dims(P,A,sigma,rho)
        K    = _kktutils_make_kkt(P,A,sigma,rho,:F)
        fact = ldlt(K)

        return new(fact,K,m,n)

    end
end

ldiv!(s::CholmodKKTSolver, rhs::Vector) = rhs .= s.fact\rhs

restart!(s::CholmodKKTSolver) = nothing #direct method, nothing to restart

function update_rho!(s::CholmodKKTSolver, rho::Union{Vector,AbstractFloat})

    didx = diagind(s.K, 0)
    @views s.K[didx[(s.n+1):(s.m+s.n)]] .= -1. ./ rho
    #complete restart
    s.fact = ldlt(s.K)
end


# -------------------------------------
# Pardiso
# -------------------------------------

abstract type AbstractPardisoKKTSolver end

#---------------------------
#Direct Solver Configuration
#---------------------------

mutable struct PardisoDirectKKTSolver <: AbstractPardisoKKTSolver

    ps::PardisoSolver
    K::SparseMatrixCSC
    m::Integer
    n::Integer
    work::Vector  #working memory since Pardiso doesn't solve in place

    function PardisoDirectKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC,sigma,rho)

        m, n, K, ps, work = _pardiso_common_init(P,A,sigma,rho)

        set_solver!(ps,Pardiso.DIRECT_SOLVER)
        pardisoinit(ps)

        # Analyze the matrix and compute a symbolic factorization.
        set_phase!(ps, Pardiso.ANALYSIS)
        pardiso(ps, K, work)

        # Compute the numeric factorization.
        set_phase!(ps, Pardiso.NUM_FACT)
        pardiso(ps, K, work)
        _pardiso_check_inertia(ps,m,n)

        ## set phase to solving for iterations
        set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)

        ##configure settings to solve in place
        set_iparm!(ps,6,1)

        return new(ps,K,m,n,work)
    end
end

#---------------------------
#Indirect Solver Configuration
#---------------------------

mutable struct PardisoIndirectKKTSolver <: AbstractPardisoKKTSolver

    ps::PardisoSolver
    K::SparseMatrixCSC
    m::Integer
    n::Integer
    work::Vector  #working memory since Pardiso doesn't solve in place

    function PardisoIndirectKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC,sigma,rho)

        m, n, K, ps, work = _pardiso_common_init(P,A,sigma,rho)

        set_solver!(ps,Pardiso.ITERATIVE_SOLVER)
        pardisoinit(ps)

        ## set phase to solving for iterations
        set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)

        ##configure settings to solve in place
        set_iparm!(ps,6,1)

        return new(ps,K,m,n,work)
    end
end

function _pardiso_common_init(P,A,sigma,rho)

    m,n  = _kktutils_check_dims(P,A,sigma,rho)
    K    = _kktutils_make_kkt(P,A,sigma,rho,:L)
    ps   = PardisoSolver()
    work = zeros(eltype(A),m+n)

    #set to symmetric indefinite
    set_matrixtype!(ps, Pardiso.REAL_SYM_INDEF)

    return m, n, K, ps, work
end


function ldiv!(s::AbstractPardisoKKTSolver, rhs::Vector)

    #we configure pardiso with iparm[6] = 1, so
    #it solves in place but still needs a work vector
    pardiso(s.ps, s.work, s.K, rhs)
    return
end

restart!(s::AbstractPardisoKKTSolver) = nothing

function update_rho!(s::AbstractPardisoKKTSolver, rho::Union{Vector,AbstractFloat})

    didx = diagind(s.K, 0)
    @views s.K[didx[(s.n+1):(s.m+s.n)]] .= -1. ./ rho

    #only refactor for the direct solver
    if s.ps.solver == Pardiso.DIRECT_SOLVER
        # Compute the numeric factorization again,
        #but skipping the analysis phase
        set_phase!(s.ps, Pardiso.NUM_FACT)
        pardiso(s.ps, s.K, s.work)
        _pardiso_check_inertia(s.ps,s.m,s.n)

        ## set back to solving phase
        set_phase!(s.ps, Pardiso.SOLVE_ITERATIVE_REFINE)
    end
end

function _pardiso_check_inertia(ps,m,n)

    num_pos_eigenvalues = get_iparm(ps, 22)
    num_neg_eigenvalues = get_iparm(ps, 23)

    num_neg_eigenvalues == m || throw("P matrix appears to have negative eigenvalues")
end
