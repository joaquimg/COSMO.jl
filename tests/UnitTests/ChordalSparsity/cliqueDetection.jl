
workspace()
include("../../../src/Solver.jl")
include("../../../../../PositiveDefiniteCompletion/Code/chompack_wrapper.jl")
using OSSDP, ChordalSparsity, GraphModule, chomWrap
using Base.Test


rng = MersenneTwister(3232)
nn = 1000

function getChordalMatrix(g)
  nv = length(g.adjacencyList)
  A = eye(nv)

  for iii = 1:nv
    for jjj in g.adjacencyList[iii]
      A[iii,jjj] = 1
    end
  end
  return A
end


function subset2(x,y)
    lenx = length(x)
    first = x[1]
    if lenx == 1
        return findnext(y, first, 1) != 0
    end
    leny = length(y)
    lim = length(y) - length(x) + 1
    cur = 1
    while (cur = findnext(y, first, cur)) != 0
        cur > lim && break
        beg = cur
        @inbounds for i = 2:lenx
            y[beg += 1] != x[i] && (beg = 0 ; break)
        end
        beg != 0 && return true
        cur += 1
    end
    false
end

function cliquesAreMaximal(cliques)
  jjj = 1
  const ind = collect(1:length(cliques))
  for c in cliques
    oInd = filter(i -> i!=jjj,ind)
    for o in cliques[oInd]
      if subset2(c,o)
        return false
      end
    end
    jjj +=1
  end
  return true
end



 @testset "Maximal cliques" begin
   for iii = 1:nn
    # create a random chordal graph
    dim = rand(rng,10:100)
    density = rand(rng,0.1:0.1:0.4)
    A = sprand(rng,dim,dim,density)
    A = A+A'
    g = Graph(A)

    # derive the corresponding chordal sparsity pattern in matrix form
    B = getChordalMatrix(g)

    # determine maximal cliques with chompack
    cliques_ref = chomWrap.getCliques(B)

    # determine maximal cliques with my routines
    sp = SparsityPattern(sparse(B));
    cliques = sp.cliques
    @test cliquesAreMaximal(cliques_ref)
    @test cliquesAreMaximal(cliques)
    println("$(iii)/$(nn) completed!")
  end
end


