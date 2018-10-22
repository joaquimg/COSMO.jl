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

function _kktutils_make_kkt(P,A,sigma,rho)

    R = length(rho)   == 1 ? ((1.)./rho)*I : Diagonal((1.)./rho)
    S = length(sigma) == 1 ? (sigma)*I : Diagonal(sigma)

    #compute the full KKT matrix
    K    = [P+S A'; A -R]
    return K

end


# -------------------------------------
# QDLDL solver
# -------------------------------------
mutable struct QDLDLKKTSolver <: AbstractKKTSolver

    fact::QDLDL.QDLDLFactorisation
    K
    m::Integer
    n::Integer

    function QDLDLKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC,sigma,rho)

        m,n = _kktutils_check_dims(P,A,sigma,rho)
        K   = _kktutils_make_kkt(P,A,sigma,rho)
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
        K    = _kktutils_make_kkt(P,A,sigma,rho)
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
mutable struct PardisoKKTSolver <: AbstractKKTSolver

    ps::PardisoSolver
    K::SparseMatrixCSC
    m::Integer
    n::Integer
    work::Vector  #working memory since Pardiso doesn't solve in place

    function PardisoKKTSolver(P::SparseMatrixCSC, A::SparseMatrixCSC,sigma,rho)

        m,n  = _kktutils_check_dims(P,A,sigma,rho)
        K    = tril(_kktutils_make_kkt(P,A,sigma,rho))
        ps   = PardisoSolver()
        work = Vector{eltype(A)}(undef,m+n)

        #--------
        # solve in steps following the Pardiso package symmetric
        # indefinite solver example
        #--------

        #set to symmetric indefinite
        set_matrixtype!(ps, Pardiso.REAL_SYM_INDEF)
        # Initialize the default settings with the current matrix type
        pardisoinit(ps)

        # Analyze the matrix and compute a symbolic factorization.
        set_phase!(ps, Pardiso.ANALYSIS)
        set_perm!(ps, amd(K))
        pardiso(ps, K, work)

        # Compute the numeric factorization.
        set_phase!(ps, Pardiso.NUM_FACT)
        pardiso(ps, K, work)
        _pardiso_check_inertia(ps,m,n)

        ## set phase to solving for iterations
        set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)

        return new(ps,K,m,n,work)

    end
end

function ldiv!(s::PardisoKKTSolver, rhs::Vector)

    #pardiso can't solve in place, so copy the
    #rhs in work, then solve into the rhs input
    copyto!(s.work,rhs)
    pardiso(s.ps, rhs, s.K, s.work)
    return
end


restart!(s::PardisoKKTSolver) = nothing #direct method, nothing to restart

function update_rho!(s::PardisoKKTSolver, rho::Union{Vector,AbstractFloat})

    didx = diagind(s.K, 0)
    @views s.K[didx[(s.n+1):(s.m+s.n)]] .= -1. ./ rho

    # Compute the numeric factorization again,
    #but skipping the analysis phase
    set_phase!(s.ps, Pardiso.NUM_FACT)
    pardiso(s.ps, s.K, s.work)
    _pardiso_check_inertia(s.ps,s.m,s.n)

    ## set back to solving phase
    set_phase!(s.ps, Pardiso.SOLVE_ITERATIVE_REFINE)

end

function _pardiso_check_inertia(ps,m,n)

    num_pos_eigenvalues = get_iparm(ps, 22)
    num_neg_eigenvalues = get_iparm(ps, 23)

    num_neg_eigenvalues == m || throw("P matrix appears to have negative eigenvalues")
end
