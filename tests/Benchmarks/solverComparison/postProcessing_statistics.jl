workspace()
include("../../../src/Solver.jl")
include("../solverComparison/Compare.jl")
 include("./LatexExport.jl")

# load data files
using  FileIO, OSSDP, Compare
using LatexExport

folderName = "BlockArrow-v06-Conference"
dir = "../../../test/resultDataFiles/SDP_Benchmark_Problems/"*folderName
results = []
for f in filter(x -> endswith(x, ".jld"), readdir(dir))
    f = split(f,".")[1]
    push!(results,String(f))
end
#filter!(x->x != "Combined",results)
metricsArr = Array{Compare.ProblemMetrics}(length(results))

# Step 1: loop over all problem types and the combined file and calculate important metrics
for iii=1:length(results)
  r = results[iii]
  data = load(dir*"/"*r*".jld")
  resData = data["resData"]

  contains(resData[1].problemType,"Combined") ? COMBINE_FLAG = true : COMBINE_FLAG = false
  # create data object to hold the metrics for this problem type
  pm = Compare.ProblemMetrics(resData[1].problemType,COMBINE_FLAG,resData[1].ind,length(resData))

  # loop over solver and compute metrics
  k = 1
  for s in resData
    solvedInd = find(x->x==:solved,s.status)
    # calculate mean
    meanIterAll = mean(s.iter)
    length(solvedInd) > 0 ? meanIterSolved = mean(s.iter[solvedInd]) : meanIterSolved = Inf

    numSolved = length(solvedInd)
    percSolved = numSolved/s.ind

    meanErr = mean(abs.(s.objVal - s.objTrue)[solvedInd])

    meanSolveTime = mean(s.solverTime)
    meanSetupTime = mean(s.setupTime)
    meanGraphTime = mean(s.graphTime)*1000
    meanIterTime = mean(s.iterTime)
    meanProjTime = mean(s.projTime)

    meanAvgIterTime = mean(s.iterTime./s.iter)*1000
    meanAvgProjTime = mean(s.projTime./s.iter)*1000

    meanNZ = mean(s.problemDim[:,3])
    sm = Compare.SolverMetrics(s.solverName,s.adaptionON,s.scalingON,meanIterAll,meanIterSolved,numSolved,percSolved,meanErr,meanSolveTime,meanSetupTime,meanGraphTime,meanIterTime,meanProjTime,meanAvgIterTime,meanAvgProjTime,meanNZ)
    pm.solverResults[k] = sm
    k+=1
  end
  metricsArr[iii] = pm
end

# permute metricsArr to have combined results at the end of array
cind = find(x-> x.combinedData,metricsArr)
if length(cind) > 0
  cind = cind[1]
  if cind == 1
    p = [collect(2:length(metricsArr));1]
  elseif cind > 1 && cind < length(metricsArr)
    p = [collect(1:cind-1);collect(cind+1:length(metricsArr));2]
  end
  permute!(metricsArr,p)
end
# Step 2: Print results to screen
println("-"^50)
println("Postprocessing statistics:")
for pm in metricsArr
  println("-"^50)
  println(">>> Problem Type: $(pm.problemType)")
  println("- Combined Data: $(pm.combinedData)")
  println("- Number of Problems: $(pm.numProblems)")
  println("- Number of Solvers: $(pm.numSolvers)")

    println("Solver Name:\tNum Solved:\t% Solved:\tIter(all):\tIter(sol.):\tError:\t\ttotal t.[s]:\tsetup t.[s]:\tgraph t.[ms]:\tfactor t.[ms]:\tavg iter t.[ms]:\tavg proj t.[ms]:")
  for sm in pm.solverResults

    if length(sm.solverName) < 10
      name = sm.solverName * " "^(10-length(sm.solverName))
    else
      name = sm.solverName
    end
    @printf("%s\t",name)
    @printf("%d\t\t",sm.numSolved)
    @printf("%.2f\t\t",sm.percSolved)
    @printf("%.2f\t\t",sm.meanIterAll)
    @printf("%.2f\t\t",sm.meanIterSolved)
    @printf("%.2f\t\t",sm.meanErr)
    @printf("%.3f\t\t",sm.meansolT)
    @printf("%.3f\t\t",sm.meansT)
    @printf("%.2f\t\t",sm.meangT)
    @printf("%.2f\t\t\t",sm.meanAvgIterTime)
    @printf("%.2f\n",sm.meanAvgProjTime)
  end
end
println("-"^50)

# Step 3: Print results to LaTeX table
resPath = dir*"/latex/"
!ispath(resPath) && mkdir(resPath)

createLatexTable(metricsArr,dir,pwd())



