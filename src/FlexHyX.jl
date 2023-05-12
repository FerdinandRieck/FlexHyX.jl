module FlexHyX
    using DifferentialEquations, NLsolve, Plots
    using TerminalLoggers
    using Dates
    using LinearAlgebra
    using SparseArrays
    import Pkg
    import JSON

    dir = dirname(@__DIR__)

    #-- Funktionen einfügen
    pfad = filter(contains(r".jl$"), readdir(dir*"/src/Funktionen/";join=true))
    include.(pfad)

    #-- Knoten einfügen
    pfad = filter(contains(r".jl$"), readdir(dir*"/src/Komponenten/Knoten/";join=true))
    include.(pfad)

    #-- Kanten einfügen
    pfad = filter(contains(r".jl$"), readdir(dir*"/src/Komponenten/Kanten/";join=true))
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

    #sol = solveNetzwerk(dir*"/src")
    #plotsol(sol)
end
