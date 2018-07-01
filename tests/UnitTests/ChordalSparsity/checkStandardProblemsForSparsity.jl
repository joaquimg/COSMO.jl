# Script to check the standard SDP problems for chordal sparsity:



workspace()
include("../../../src/Solver.jl")
include("../../Benchmarks/solverComparison/Compare.jl")

using OSSDP,Compare, Base.Test, JLD, OSSDPTypes, ChordalSparsity

SAVE_ALWAYS = true
maxIter = 2000


# find available problem types in SDP_Benchmark_Problems folder
probPath = "/Users/Micha/Dropbox/Research/SDP_Benchmark_Problems/DataFiles/Julia/"
existingFolders = readdir(probPath)
problemTypes = []
for f in filter(x -> !startswith(x, "."), readdir(probPath))
    f = split(f,".")[1]
    push!(problemTypes,String(f))
end
filter!(x->!in(x,["SmallestCircle"]),problemTypes)
println(">>> $(length(problemTypes)) Problem Type(s) detected!")

# run tests for each problem type
for pType in problemTypes
  println(">>> ProblemType: $(pType)")
  dirPath = probPath*pType*"/"

  # find all file names in problem type folder
  problems = []
  for f in filter(x -> endswith(x, ".jld"), readdir(dirPath))
      f = split(f,".")[1]
      push!(problems,String(f))
  end
  nn = length(problems)


  # loop over all problems in problem type folder
  for iii =1:1:length(problems)
    iii == length(problems) && (SAVE_ALWAYS = true)
    gc()
    problem = problems[iii]
    data = JLD.load("$(dirPath)"*"$(problem).jld")
    P = sparse(data["P"])
    q = full(data["q"])[:]
    r = data["r"]
    A = sparse(data["A"])
    b = full(data["b"])[:]
    m = data["m"]
    n = data["n"]
    Kf = data["Kf"]
    Kl = data["Kl"]
    Kq = data["Kq"]
    Ks = data["Ks"]

    # define cone membership
    K = Cone(Kf,Kl,Kq,Ks)
    settings = OSSDPSettings(rho=0.1,sigma=1e-6,alpha=1.6,max_iter=3,verbose=false,adaptive_rho=true,scaling=10,decompose=true)

    chordalInfo = OSSDPTypes.ChordalInfo(size(A,1),size(A,2),K)
    Pa, qa, Aa, ba, Ka = ChordalSparsity.chordalDecomposition!(P,q,A,b,K,settings,chordalInfo)
    if length(Ka.s) > length(K.s)
      println("Problem $(problem): Number of cones increased from $(length(K.s)) to $(length(Ka.s)).")
    end
  end
  println(">>> ProblemType: $(pType) completed!")
end
