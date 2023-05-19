function plotsol(sol,x,y)
    p1 = plot(sol.t/3600,sol'[:,x], linewidth = 2, xlabel = "Zeit /h") # ylims=(0.494, 0.5)
    p2 = plot(sol.t/3600,sol'[:,y]/65880, linewidth = 2, xlabel = "Zeit /h", ylabel = "soc", title = "Ladezustand", label = "soc",legend=:topleft)
    p = plot(p1,p2,layout = (2, 1))
    display(p)
end