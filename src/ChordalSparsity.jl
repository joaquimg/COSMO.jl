module ChordalSparsity


using GraphModule, TreeModule
export SparsityPattern

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


end #MODULE