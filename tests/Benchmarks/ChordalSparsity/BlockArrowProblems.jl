# Test routine that solve all SDP Benchmark problems and saves the results in .jld resultDataFiles
# uses the Compare Module

workspace()
include("../../../src/Solver.jl")
include("../solverComparison/Compare.jl")

using OSSDP, Base.Test, Compare, JLD
SAVE_ALWAYS = true
maxIter = 2000

# choose experiment name, otherwise save in folder based on timestamp
EXPERIMENT_NAME = "BlockArrowProblems"
# create a new folder for the results based on timestamp
timestamp = Dates.format(now(), "yyddmm_HH-MM")

isdefined(:EXPERIMENT_NAME) ? folderName=EXPERIMENT_NAME : folderName=timestamp
resPath = "../resultDataFiles/SDP_Benchmark_Problems/"*folderName
!ispath(resPath) && mkdir(resPath)

# find available problem types in SDP_Benchmark_Problems folder
# probPath = "/Users/Micha/Dropbox/Research/SDP_Benchmark_Problems/DataFiles/Julia/"
probPath = "/Users/Micha/Dropbox/Research/SDP_Benchmark_Problems/DataFiles/BlockArrow/"
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
  resFileName = resPath*"/"*pType*".jld"
  dirPath = probPath*pType*"/"

  # find all file names in problem type folder
  problems = []
  for f in filter(x -> endswith(x, ".jld"), readdir(dirPath))
      f = split(f,".")[1]
      push!(problems,String(f))
  end
  nn = length(problems)


  problemType = JLD.load("$(dirPath)"*"$(problems[1]).jld","problemType")
  sr1 = SolverResult(nn, problemType,"QOCS - normal",timestamp,0,true,true)
  sr2 = SolverResult(nn, problemType,"QOCS - decomp",timestamp,0,true,true)
  sr3 = SolverResult(nn, problemType,"MOSEK",timestamp,0,true,false)

  resData = [sr1;sr2;sr3]

  # loop over all problems in problem type folder
  for iii =1:1:length(problems)
    iii == length(problems) && (SAVE_ALWAYS = true)
    gc()
    problem = problems[iii]
    data = JLD.load("$(dirPath)"*"$(problem).jld")
    P = data["P"]
    q = data["q"]
    r = data["r"]
    A = data["A"]
    b = data["b"]
    m = data["m"]
    n = data["n"]
    Kf = data["Kf"]
    Kl = data["Kl"]
    Kq = data["Kq"]
    Ks = data["Ks"]

    # objTrue = data["objTrue"]
    problemName = data["problemName"]


    pDims = [size(A,1);size(A,2);nnz(A)]

    settings            = OSSDPSettings(max_iter=maxIter,checkTermination=40,adaptive_rho=true,decompose=false)
    settings_decompose  = OSSDPSettings(max_iter=maxIter,checkTermination=40,adaptive_rho=true,decompose=true)

    # define cone membership
    K = Cone(Kf,Kl,Kq,Ks)

    # Solve with QOCS
    res1,nothing = OSSDP.solve(P,q,A,b,K,settings);
    print("\n.")
    res2,nothing = OSSDP.solve(P,q,A,b,K,settings_decompose);
    print("\n.")

    # solve with MOSEK


    # solve with SCS v.1.2.6


    # save results
    # resArray = [res1;res2;res3;res4;res5;res6;res7;res8]
    resArray = [res5;res6;res7;res8]
    updateResults!(resFileName,resData,resArray,pDims,problemName,r,SAVE_ALWAYS)
    printStatus(iii,nn,problemName,resData)
  end
  println(">>> ProblemType: $(pType) completed!")
end

