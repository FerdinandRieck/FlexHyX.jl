module FlexHyX
    using DifferentialEquations, NLsolve, Plots
    using TerminalLoggers
    using Dates
    using LinearAlgebra
    using SparseArrays
    import JSON

    dir = dirname(@__FILE__)

    #-- Funktionen einfügen
    pfad = filter(contains(r".jl$"), readdir(dir*"/Funktionen/";join=true))
    include.(pfad)

    #-- Knoten einfügen
    pfad = filter(contains(r".jl$"), readdir(dir*"/Komponenten/Knoten/";join=true))
    include.(pfad)

    #-- Kanten einfügen
    pfad = filter(contains(r".jl$"), readdir(dir*"/Komponenten/Kanten/";join=true))
    include.(pfad)

    #=
    #-- Typenhierarchie anzeigen
    using AbstractTrees
	AbstractTrees.children(x::Type) = subtypes(x)
    print_tree(FlexHyX)
    =#

    #-- Funktionen exportieren
    export solveNetzwerk
    export plotsol

    idx_ele, sol, y = solveNetzwerk(dir)

    #=
    iV3 = idx_ele["3i"][2] 
    iPV4 = idx_ele["4i"][2]
    iB1 = idx_ele["1i"][2]
    iS7 = idx_ele["7i"][2]
    iSP8 = idx_ele["8i"][2]
    iD18 = idx_ele["18i"][2]
    iBZ12 = idx_ele["12i"][2]
    iE9 = idx_ele["9i"][2]
    U3 = idx_ele["14P"][2]
    SOC = idx_ele["1q"][2]
    Θ = idx_ele["15Θ"][2]
    plotsol(sol,iS7,SOC)
    plotsol(y,sol.t,48,28)
    plotsol(sol)
    =#
    plotsol(y,sol.t)
 end

#=
Aktuelle Bugs:
- Eventfunktion muss aktuell noch in solveNetzwerk Funktion oben includiert werden
- leichte unterschiede der Eventzeitpunkte von Matlab und JULIA mit Zeitreihe
=#