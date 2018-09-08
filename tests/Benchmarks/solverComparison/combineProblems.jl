# workspace()
# include("../../../src/Solver.jl")
# include("../solverComparison/Compare.jl")

# using  Compare, OSSDP, FileIO
folderName = "BlockArrow-v06-Conference"
dir = "../../../test/resultDataFiles/SDP_Benchmark_Problems/"*folderName
results = []
for f in filter(x -> endswith(x, ".jld"), readdir(dir))
    f = split(f,".")[1]
    push!(results,String(f))
end
filter!(x->!in(x,["Combined"]),results)

length(unique(results)) != length(results) && warn("Careful, you are trying to combine files from the same problem type.")


# determine how many solvers are compared and determine number of total problems
problemSizes = zeros(length(results))
numSolvers = 0
for iii = 1:length(results)
  data = load(dir*"/"*results[iii]*".jld")
  rData = data["resData"]
  problemSizes[iii] = rData[1].ind
  iii == 1 && (numSolvers = length(data["resData"]))
end
numProblems = Int(sum(problemSizes))

# create array that holds combined data, each entry represents one solver
resCombined = Array{Compare.SolverResult}(numSolvers)


# MAIN LOOP OVER PROBLEM TYPES
for iii=1:length(results)

  r = results[iii]
  data = load(dir*"/"*r*".jld")
  resData = data["resData"]
  length(resData) != numSolvers && error("Number of solvers in $(r) doesnt match other files ( $(length(resData)) vs. $(numSolvers) ).")

  # for the first problem type recreate the solver result objects
  if iii == 1
    # create new solver results object that can then be filled with combined data
    for k = 1:numSolvers
      problemDim_size = size(resData[k].problemDim,2)
      resCombined[k] = SolverResult(numProblems, "Combined problems",resData[k].solverName,resData[k].timeStamp,0,resData[k].scalingON,resData[k].adaptionON)
    end
  end

  # append results of individual problems to combined results
  for k=1:length(resData)
    resSource = resData[k]
    resTarget = resCombined[k]
    if resSource.solverName != resTarget.solverName || resSource.adaptionON != resTarget.adaptionON || resSource.scalingON != resTarget.scalingON
      error("You are trying to combine the data of two different solver settings: $(resTarget.solverName) vs. $(resSource.solverName).")
    end
    # add check that solver is still the same

    resTarget.iter[resTarget.ind+1:resTarget.ind+resSource.ind] = resSource.iter
    resTarget.objVal[resTarget.ind+1:resTarget.ind+resSource.ind] = resSource.objVal
    resTarget.solverTime[resTarget.ind+1:resTarget.ind+resSource.ind] = resSource.solverTime
    resTarget.setupTime[resTarget.ind+1:resTarget.ind+resSource.ind] = resSource.setupTime
    resTarget.graphTime[resTarget.ind+1:resTarget.ind+resSource.ind] = resSource.graphTime
    resTarget.iterTime[resTarget.ind+1:resTarget.ind+resSource.ind] = resSource.iterTime
    resTarget.projTime[resTarget.ind+1:resTarget.ind+resSource.ind] = resSource.projTime
    resTarget.status[resTarget.ind+1:resTarget.ind+resSource.ind] = resSource.status
    resTarget.problemDim[resTarget.ind+1:resTarget.ind+resSource.ind,:] = resSource.problemDim
    resTarget.problemName[resTarget.ind+1:resTarget.ind+resSource.ind] = resSource.problemName
    resTarget.ind += resSource.ind
  end


end
fn = dir * "/Combined.jld"
save(fn, "resData", resCombined)

