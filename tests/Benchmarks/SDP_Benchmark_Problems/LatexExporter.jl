module LatexExport

using Compare
export createLatexTable, createIterLatexTable


function createLatexTable(metricsArr::Array{Compare.ProblemMetrics},dir::String)


  caption1 = "Performance of different solver configurations over the set of SDP Benchmark problems (mean values)."

  open(dir*"/table.tex", "w") do file
    # create the header of the table
    write(file, "\\begin{longtable}{l c c c c c c c c c c}\n\\toprule\n")
    # add data
    for pm in metricsArr
      write(file, "\\multicolumn{7}{l}{Problem Type: \\textbf{$(pm.problemType)}} \\\\ \n \\multicolumn{7}{l}{Number of Problems: $(pm.numProblems)} \\\\ \n \\multicolumn{7}{l}{Number of Solvers: $(pm.numSolvers)} \\\\ \n  \\hline\n")
      # write(file, " Solver Name &  \\# solved & \\% solved & NZ & Iter(solv) & RT (all) & IT (all) & avgIterT (all) & setupT (all)  & Error(solv) \\\\ [0.5ex]\n\\hline\n")
      write(file, " Solver Name &  \\# solved & \\% solved & NZ & Iter(solv) & RT (all) & IT (all) & avgIterT (all) & setupT (all)   \\\\ [0.5ex]\n\\hline\n")
      for sm in pm.solverResults
        # write(file, "$(sm.solverName)  & $(sm.numSolved) & $(round(sm.percSolved,2)) & $(round(sm.meanNZ,2))  & $(round(sm.meanIterSolved,2)) & $(round(sm.meanRunTime,4)) & $(round(sm.meanIterTime,4))  & $(round(sm.meanAvgIterTime,4))& $(round(sm.meanSetupTime,4)) & $(round(sm.meanErr,4))\\\\ \n")
        write(file, "$(sm.solverName)  & $(sm.numSolved) & $(round(sm.percSolved,2)) & $(round(sm.meanNZ,2))  & $(round(sm.meanIterSolved,2)) & $(round(sm.meanRunTime,4)) & $(round(sm.meanIterTime,4))  & $(round(sm.meanAvgIterTime,4))& $(round(sm.meanSetupTime,4)) \\\\ \n")
      end
      write(file,"\\bottomrule \\\\ \n")
    end
    # write footer
    write(file, "\\bottomrule \n\\caption{$(caption1)}\n\\label{table:1}\n\\end{longtable}")

  end
  return nothing
end

  function createIterLatexTable(resComp,numProblemsArr,results,dir)
  caption1 = "Number of problems where a solver configuration with scaling performed better than the unscaled configuration (only cases with rho-adaption)."
  open(dir*"/iterTable.tex", "w") do file
    write(file, "\\begin{longtable}{l c c c }\n\\toprule\n")
    for iii=1:length(numProblemsArr)
      numProblems = numProblemsArr[iii]
      write(file, "\\multicolumn{4}{l}{Problem Type: \\textbf{$(results[iii])}} \\\\ \n \\multicolumn{4}{l}{Number of Problems: $(numProblemsArr[iii])} \\\\ \n  \\hline\n")
      write(file, " Mean \$ \\leq \$ Unscaled &  Geom \$ \\leq \$ Unscaled & Sym \$ \\leq \$ Unscaled &  Any Scaled \$ \\leq \$  Unscaled \\\\ [0.5ex]\n")
      write(file, "$(resComp[iii,1]) ($(round(resComp[iii,1]/numProblems*100,2)) \\%) & $(resComp[iii,2]) ($(round(resComp[iii,2]/numProblems*100,2)) \\%) & $(resComp[iii,3]) ($(round(resComp[iii,3]/numProblems*100,2)) \\%)  & $(resComp[iii,4]) ($(round(resComp[iii,4]/numProblems*100,2)) \\%)  \\\\ \n")
      write(file,"\\bottomrule \\\\ \n")
    end
        write(file, " \n\\caption{$(caption1)}\n\\label{table:2}\n\\end{longtable}")

  end
  return nothing

  end


 # export settings object as latex table
  function exportSolverSettings(settings,dir)
  caption1 = "Solver Settings"
  fn = fieldnames(settings)
  sort!(fn)
  cols = 4
  n = length(fn)
  rows = floor(n/cols)
  leftovers = n % cols
  ind = collect(1:n)


  open(dir*"/solverSettings.tex", "w") do file
    write(file, "\\begin{longtable}{$("r c "^cols)}\n\\toprule\n")
    while !isempty(ind)
      iii = ind[1]
      name = replace(string(fn[iii]), "_","\\_")

        if mod(iii,cols) == 1
          write(file, " $(name)  &  $(getfield(settings,fn[iii])) " )
        else
          write(file, " & $(name)  &  $(getfield(settings,fn[iii])) " )
        end
        (mod(iii,cols) == 0) && write(file,"\\\\ [0.5ex]\n")
        deleteat!(ind,1)
    end
    write(file, "\\\\ \n \\bottomrule \n\\caption{$(caption1)}\n\\label{table:2}\n\\end{longtable}")

  end
  return nothing

  end
end #MODULE