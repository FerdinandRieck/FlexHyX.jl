function plotSol(sol,x,y)
    p = plot(sol.t,sol'[:,x], linewidth = 2, xlabel = "Zeit /Tage", xlims=(0,40000))
    #p = plot(sol.t/(3600*24),sol'[:,x], linewidth = 2, xlabel = "Zeit /Tage")#, xlims=(800/(3600*24), 1030/(3600*24)))
    #p2 = plot(sol.t/3600,sol'[:,y]/(18.3*3600*n_par), linewidth = 2, xlabel = "Zeit /h", ylabel = "soc", title = "Ladezustand", label = "soc",legend=:topleft)
    #p = plot(p1,p2,layout = (2, 1))
    #p = plot(sol,idxs=x, linewidth = 2, xlabel = "Zeit /s")#, ylims=(101300, 103000))
    return p
end

function plotSol(sol,t,x,y)
  n_par = 1
  p1 = plot(t/3600,sol[x,:], linewidth = 2, xlabel = "Zeit /h") # ylims=(0.494, 0.5)
  p2 = plot(t/3600,sol[y,:]/(18.3*3600*n_par), linewidth = 2, xlabel = "Zeit /h", ylabel = "soc", title = "Ladezustand", label = "soc",legend=:topleft)
  p = plot(p1,p2,layout = (2, 1))
  return p
end

function plotSol(sol)
  n_par = 1
  p1 = plot(sol.t/3600,sol'[:,[29,30]], linewidth = 2, xlabel = "Zeit /h", ylabel = "A", title = "Ströme", label = ["iV" "iPV"],legend=:bottomleft, legend_columns=-1)
  p2 = plot(sol.t/3600,sol'[:,[26,34,39]], linewidth = 2, xlabel = "Zeit /h", ylabel = "A", title = "Ströme", label = ["iB" "iE" "iBZ"],legend=:bottom, legend_columns=-1)
  p3 = plot(sol.t/3600,sol'[:,[2]], linewidth = 2, xlabel = "Zeit /h", ylabel = "V", title = "Spannungen", label = "U3", legend=:bottom)
  p4 = plot(sol.t/3600,sol'[:,28]/(18.3*3600*n_par), linewidth = 2, xlabel = "Zeit /h", ylabel = "%", title = "Ladezustand", label = "soc",legend=:bottom)
  p5 = plot(sol.t/3600,sol'[:,22], linewidth = 2, xlabel = "Zeit /h", ylabel = "%", title = "Beladung MHS", label = "Θ",legend=:top)
   
  p = plot(p1,p2,p3,p4,p5, layout = (5, 1), size=(1000,1000))
  return p
end

function plotSol(sol,t)
  n_par = 1
  p1 = plot(t/3600,sol[50,:], linewidth = 2, xlabel = "Zeit /h", ylabel = "W", title = "Leistung", label = "3iV-P",legend=:topright, legend_columns=-1)
  p1 = plot!(t/3600,sol[51,:], linewidth = 2, label = "4iPV-P")
  p2 = plot(t/3600,sol[26,:], linewidth = 2, xlabel = "Zeit /h", ylabel = "A", title = "Ströme", label = "1iB",legend=:bottom, legend_columns=-1)
  p2 = plot!(t/3600,sol[34,:], linewidth = 2, label = "9iE")
  p2 = plot!(t/3600,sol[39,:], linewidth = 2, label = "iBZ")
  p3 = plot(t/3600,sol[2,:], linewidth = 2, xlabel = "Zeit /h", ylabel = "V", title = "Spannung", label = "U3", legend=:bottom)
  p4 = plot(t/3600,sol[28,:]/(18.3*3600*n_par), linewidth = 2, xlabel = "Zeit /h", ylabel = "%", title = "Ladezustand", label = "soc",legend=:bottom)
  p5 = plot(t/3600,sol[22,:], linewidth = 2, xlabel = "Zeit /h", ylabel = "%", title = "Beladung MHS", label = "Θ",legend=:top)
  
  p = plot(p1,p2,p4,p5, layout = (4, 1), size=(1000,1000))
  return p
end