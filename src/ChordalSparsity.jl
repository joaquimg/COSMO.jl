module ChordalSparsity


using GraphModule, TreeModule, OSSDPTypes
export SparsityPattern, findStackingMatrix, CliqueSet, getClique, chordalDecomposition!, reverseDecomposition!
export vecToMatInd,findCommonSparsityPattern, findStackingMatrix


# ---------------------------
# STRUCTS
# ---------------------------
mutable struct SparsityPattern
  g::Graph
  cliqueTree::Tree
  cliques::Array{Array{Int64,1}}


  # constructor for sparsity pattern
  function SparsityPattern(A::SparseMatrixCSC{Int64,Int64})
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
A1 = [1. 2 0; 4 5 6; 0 8 9]
A2 = [1. 2 0; 4 5 6; 0 8 9]
A = sparse([[9; 9; 9];vec(A1);vec(A2)])[:,:]
b = rand(21,1)
b[[6,10,15,19]] = 0
K = OSSDPTypes.Cone(1,2,[],[9 9])
P = speye(size(A,2))
q = zeros(size(A,2))

function chordalDecomposition!(P::SparseMatrixCSC{Float64,Int64},q::Vector{Float64},A::SparseMatrixCSC{Float64,Int64},b::Vector{Float64},K::OSSDPTypes.Cone,settings::OSSDPTypes.OSSDPSettings,chordalInfo::OSSDPTypes.ChordalInfo)

  # do nothing if no psd cones present in the problem
  if length(K.s) == 0
    settings.decompose = 0
    return P,q,A,b,K
  end

  # find sparsity pattern for each cone
  numCones = length(K.s)
  numCones ==0 && return nothing
  spArr = Array{ChordalSparsity.SparsityPattern,1}(numCones)
  cliqueSets = Array{ChordalSparsity.CliqueSet,1}(numCones)
  st = K.f+K.l+sum(K.q) + 1

  # find sparsity pattern graphs and clique sets for each cone
  for iii=1:numCones
    e =st+K.s[iii] - 1
    csp = findCommonSparsityPattern(A[st:e,:],b[st:e])
    sp = SparsityPattern(csp)
    spArr[iii] = sp
    cliqueSets[iii] = CliqueSet(sp.cliques)
    st+= K.s[iii]
  end

  # find transformation matrix H and store it
  H, KsA = findStackingMatrix(K,cliqueSets)
  chordalInfo.H = H

  # augment the system, change P,q,A,b
  m,n = size(A)
  mH,nH = size(H)
  P = blkdiag(P,spzeros(nH,nH))
  q = [q;zeros(nH)][:]
  A = [A H; spzeros(nH,n) -speye(nH,nH)]
  b = [b;zeros(nH)][:]

  K = OSSDPTypes.Cone(K.f + mH,K.l,K.q,KsA)
  # problem data is returned instead of mutated in the function since problem size changes
  # TODO: Also return some information about sparsity to reverse decomposition
  return P,q,A,b,K
end


# Asub = sparse(rand(16,4))
# bsub = rand(20,1)

# Asub[4,:] = 0
# Asub[15,:] = 0
# Asub[8,:] = 0
# bsub[15] = 0
# bsub[8] = 0

function findCommonSparsityPattern(Asub,bsub)
  m,n = size(Asub)
  AInd = find(x->Asub[x,:] == zeros(n),collect(1:m))
  bInd = find(x->x==0,bsub)
  zeroRows = intersect(AInd,bInd)
  mSize = Int(sqrt(m))
  csp = spzeros(Int64,mSize,mSize)
  csp[:,:] = 1

  for ind in zeroRows
    i,j = vecToMatInd(ind,mSize)
    csp[i,j] = 0
  end
  return csp
end

function vecToMatInd(ind::Int64,n::Int64)
  ind > n^2 && error("Index ind out of range.")
  ind == 1 && (return 1,1)

  r = ind % n

  if r == 0
    j = Int(ind/n)
    i = n
  else
    j = Int(floor(ind/n) + 1)
    i = r
  end
  return i,j
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
  n = bK + sum(stackedSizes)
  # length of original vector s
  m = bK + sum(K.s)

  H = spzeros(m,n)
  H[1:bK,1:bK] = speye(bK)
  bK += 1
  b = bK
  Ks = Array{Int64}(0)
  for kkk = 1:length(cliqueSets)
    cliqueSet = cliqueSets[kkk]
    mH = Int64
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
      Hkt = (kron(Ek,Ek))'
      mH,nH = size(Hkt)
      H[b:b+mH-1,bK:bK+nH-1] = Hkt
      bK += nH
      push!(Ks,nH)
    end
    b += mH
  end
  return H, Ks

end

function reverseDecomposition!(ws,settings)
  mO = ws.ci.originalM
  nO = ws.ci.originalN
  H = ws.ci.H
  ws.s = H*ws.s[mO+1:end]
  ws.x = ws.x[1:nO]
  ws.μ = H*ws.μ[mO+1:end]
  # positive semidefinite completion on parts of μ
  # settings.completeDual && psdCompletion(ws.μ)
  ws.ν = H*ws.ν[mO+1:end]

  return H
end




end #MODULE