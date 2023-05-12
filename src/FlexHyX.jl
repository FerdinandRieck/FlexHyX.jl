module FlexHyX
using Plots, NLsolve, DifferentialEquations

    import Pkg
    import JSON
    using TerminalLoggers
    using Plots
    using DifferentialEquations
    using Dates
    using NLsolve
    using LinearAlgebra
    using SparseArrays

    #-- Funktionen einfügen
    dir = filter(contains(r".jl$"), readdir("src/Funktionen/";join=true))
    dir = chop.(dir, tail=0,head = 4)
    include.(dir)

    #-- Knoten einfügen
    dir = filter(contains(r".jl$"), readdir("src/Komponenten/Knoten/";join=true))
    dir = chop.(dir, tail=0,head = 4)
    include.(dir)

    #-- Kanten einfügen
    dir = filter(contains(r".jl$"), readdir("src/Komponenten/Kanten/";join=true))
    dir = chop.(dir, tail=0,head = 4)
    include.(dir)

    #=
    #-- Typenhierarchie anzeigen
    using AbstractTrees
	AbstractTrees.children(x::Type) = subtypes(x)
    print_tree(FlexHyX)
    =#

    #-- Funktionen exportieren
    export solveNetzwerk
    export plotsol

    sol = solveNetzwerk()
    plotsol(sol)
end

#Pkg.add(Pkg.PackageSpec(;name="DifferentialEquations", version="7.5.0"))