# workspace()
# include("../../../src/Solver.jl")
# include("../solverComparison/Compare.jl")

# # load data files
# using JLD, OSSDP, Compare, PyPlot

cc =  ["#1f77b4" "#ff7f0e" "#2ca02c" "#d62728" "#9467bd" "#8c564b"]
markerList = ["o" "^" "v" "s" "+" "d" "*" "8"]
xlabels = ["d - block size";"l - number of blocks";"m - number of constr" ]
timestamp = Dates.format(now(), "yyddmm_HH-MM")
  ind = [6;4;5]


folderName = "BlockArrow-v06-Conference-all"
dir = "../../../test/resultDataFiles/SDP_Benchmark_Problems/"*folderName
results = []
for f in filter(x -> endswith(x, ".jld"), readdir(dir))
    f = split(f,".")[1]
    push!(results,String(f))
end
filter!(x->!in(x,["Combined"]),results)

for iii=1:length(results)
  r = results[iii]
  resData = JLD.load(dir*"/"*r*".jld")["resData"]


  PyPlot.figure(iii,facecolor="white",figsize=(15,7))
  kkk = 1
  xdata = 1
  for s in resData
    xdata=s.problemDim[1:s.ind,ind[iii]]
    label = replace(s.solverName,"QOCS","OSCP")
    label = replace(label,"- normal","")
    PyPlot.semilogy(xdata,s.solverTime[1:s.ind],color=cc[kkk],label=label,marker=markerList[kkk])
    kkk+=1
  end
  PyPlot.grid(true)
  PyPlot.xticks(xdata,fontsize=15)
  PyPlot.yticks(fontsize=15)
  PyPlot.title(resData[1].problemType, loc="left")
  PyPlot.xlabel(xlabels[iii],fontsize=15)
  PyPlot.ylabel("Runtime [s]",fontsize=15)
  PyPlot.legend(ncol = 2,bbox_to_anchor=(0., 1.05, 1., .102),fontsize=15)
  PyPlot.savefig(dir*"/latex/$(timestamp)-a.eps")

end