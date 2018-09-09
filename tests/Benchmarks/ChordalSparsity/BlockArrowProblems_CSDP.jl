# Test routine that solve all SDP Benchmark problems and saves the results in .jld resultDataFiles
# uses the Compare Module

workspace()
include("../../../src/Solver.jl")
include("../solverComparison/Compare.jl")
include("../solverComparison/ExternalSolver.jl")

using Base.Test, Compare, JLD, ExternalSolver, CSDP, OSSDP, JuMP
SAVE_ALWAYS = true
maxIter = 2000

# choose experiment name, otherwise save in folder based on timestamp
EXPERIMENT_NAME = "BlockArrow-v06-Conference-CSDP/"
# create a new folder for the results based on timestamp
timestamp = Dates.format(now(), "yyddmm_HH-MM")

isdefined(:EXPERIMENT_NAME) ? folderName=EXPERIMENT_NAME : folderName=timestamp
resPath = "/Users/Micha/Dropbox/Research/OSSDP/Code/test/resultDataFiles/SDP_Benchmark_Problems/"*folderName
!ispath(resPath) && mkdir(resPath)

# find available problem types in SDP_Benchmark_Problems folder
# probPath = "/Users/Micha/Dropbox/Research/SDP_Benchmark_Problems/DataFiles/Julia/"
probPath = "/Users/Micha/Dropbox/Research/SDP_Benchmark_Problems/DataFiles/DecomposableProblems/BlockArrow/"
existingFolders = readdir(probPath)
problemTypes = []
# for f in filter(x -> !startswith(x, "."), readdir(probPath))
#     f = split(f,".")[1]
#     push!(problemTypes,String(f))
# end
# filter!(x->!in(x,["varD";"varM"]),problemTypes)
problemTypes = ["varL";"varD"]
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
  sr1 = SolverResult(nn, problemType,"CSDP",timestamp,0,true,true)


  resData = [sr1]

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

    objTrue = data["objTrue"]
    solveTime = data["solveTime"]
    problemName = data["problemName"]


     mconstr = data["mconstr"]
    l = data["l"]
    d = data["d"]

    pDims = [size(A,1);size(A,2);nnz(A);l;mconstr;d]

    # solve with CSDP
    m,n = size(A)
    cDim = Int(sqrt(m))
    model = Model(solver=CSDPSolver())
    @variable(model, x[1:n])
    @variable(model, S[1:cDim,1:cDim],SDP)
    s = vec(S)
    @objective(model, Min,q'*x)
    @constraint(model, A*x.+s .== b)
    solveTime = @elapsed status = JuMP.solve(model)
    objVal = getobjectivevalue(model)
    res1 = OSSDPTypes.OSSDPResult(zeros(n),zeros(m),zeros(m),zeros(m),objVal,1,:Solved,solveTime,0.,0.,0.,0.,0.,0.);

    # save results
    # resArray = [res1;res2;res3;res4;res5;res6;res7;res8]
    resArray = [res1]
    updateResults!(resFileName,resData,resArray,pDims,problemName,r,SAVE_ALWAYS,objTrue)
    printStatus(iii,nn,problemName,resData)
  end
  println(">>> ProblemType: $(pType) completed!")
end

