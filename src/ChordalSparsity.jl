module ChordalSparsity


using GraphModule, TreeModule, OSSDPTypes
export SparsityPattern, findStackingMatrix, CliqueSet, getClique, chordalDecomposition!, reverseDecomposition!



# ---------------------------
# STRUCTS
# ---------------------------
mutable struct SparsityPattern
  g::Graph
  cliqueTree::Tree
  cliques::Array{Array{Int64,1}}


  # constructor for sparsity pattern
  function SparsityPattern(A::SparseMatrixCSC{Float64,Int64})
    g = Graph(A)
    elimTree = createTreeFromGraph(g)
    superNodeElimTree = createSupernodeEliminationTree(elimTree,g)
    cliqueTree = createCliqueTree(superNodeElimTree,g)


    # for now: take cliques from cliqueTree
    cliques = Array{Array{Int64,1}}(length(cliqueTree.nodes))
    iii = 1
    for node in cliqueTree.nodes
      cliques[iii] = sort(union(node.value_top,node.value_btm))
      iii +=1
    end

    return new(g,cliqueTree,cliques)
  end

end



mutable struct CliqueSet
  cliqueInd::Vector{Int64} # an index set pointing to first element of each clique
  cliques::Vector{Int64} # just a vector with all cliques stacked
  vlen::Int64
  nBlk::Vector{Int64} # sizes of blocks
  N::Int64 #number of cliques in the set

  function CliqueSet(cliqueArr::Array{Array{Int64,1}})
    p = length(cliqueArr)
    # initialize fields
    cliqueInd = zeros(p)
    cliques = zeros(Int64,sum(map(x->length(x),cliqueArr)))
    nBlk = zeros(p)
    b = 1
    iii = 1
    for c in cliqueArr
      cliqueInd[iii] = b
      cliques[b:b+length(c)-1] = c
      nBlk[iii] = length(c)^2
      b +=length(c)
      iii +=1
    end

    return new(cliqueInd,cliques,length(cliques),nBlk,p)
  end
end



# ---------------------------
# FUNCTIONS
# ---------------------------

function chordalDecomposition!(A,b,K)

  # find sparsity pattern for each cone


  # find graphs and clique sets for each cone


  # find transformation matrix H

  # augment the system, change P,q,A,b


end

function reverseDecomposition!(ws)
  return nothing
end

function getClique(cs::ChordalSparsity.CliqueSet,ind::Int64)
  len = length(cs.cliqueInd)
  ind > len && error("Clique index ind=$(ind) is higher than number of cliques in the provided set:$(len).")
  ind < len ? (c = cs.cliques[cs.cliqueInd[ind]:cs.cliqueInd[ind+1]-1]) : (c = cs.cliques[cs.cliqueInd[ind]:end])
  return c
end

# function finds the transformation matrix H to decompose the vector s into its parts and stacks them into sbar, also returns the decomposed Ks
function findStackingMatrix(K::OSSDPTypes.Cone,cliqueSets::Array{ChordalSparsity.CliqueSet,1})

  length(K.s) != length(cliqueSets) && error("Length of K.s and number of clique sets don't match.")

  numCones = length(K.s)
  stackedSizes = zeros(Int64,numCones)
  for iii=1:numCones
    stackedSizes[iii] = sum(cliqueSets[iii].nBlk)
  end

  bK = K.f+K.l+sum(K.q)
  # length of stacked vector sBar
  m = bK + sum(stackedSizes)
  # length of original vector s
  n = bK + sum(K.s)

  H = spzeros(m,n)
  H[1:bK,1:bK] = speye(bK)
  bK += 1
  b = bK
  Ks = Array{Int64}(0)
  for kkk = 1:length(cliqueSets)
    cliqueSet = cliqueSets[kkk]
    nH = Int64
    for iii=1:cliqueSet.N
      mk = Int(sqrt(cliqueSet.nBlk[iii]))
      nk = Int(sqrt(K.s[kkk]))
      Ek = zeros(mk,nk)
      c = getClique(cliqueSet,iii)
      jjj = 1
      for v in c
        Ek[jjj,v] = 1
        jjj+=1
      end
      Hk = kron(Ek,Ek)
      mH,nH = size(Hk)
      H[b:b+mH-1,bK:bK+nH-1] = Hk
      b += mH
      push!(Ks,mH)
    end
    bK += nH
  end
  return H, Ks

end

end #MODULE