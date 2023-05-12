function plotsol(sol)
  #  p1 = plot(sol.t/3600,sol'[:,[1]], linewidth = 2, xlabel = "Zeit /h", ylabel = "V", title = "Spannungen", label = ["U0" "U1" "uB"],legend=:inside)
    #p2 = plot(sol.t/3600,sol'[:,[33]], linewidth = 2, xlabel = "Zeit /h", ylabel = "A", title = "Ströme", label = ["iV" "iPV" "iB"],legend=:bottom)
    p3 = plot(sol.t/3600,sol'[:,30]/65880, linewidth = 2, xlabel = "Zeit /h", ylabel = "soc", title = "Ladezustand", label = "soc",legend=:topleft)
   # p4 = plot(sol.t/3600,sol'[:,[17]], linewidth = 2, xlabel = "Zeit /h", ylabel = "A", title = "Öffnungsgrad", label = ["A"],legend=:inside);

    #p1 = plot(sol.t/3600,sol'[:,[7]]/65880, linewidth = 2, xlabel = "Zeit /h", ylabel = "SOC", title = "SOC", label = ["SOC"],legend=:inside)
    p = plot(p3)
    display(p)
end