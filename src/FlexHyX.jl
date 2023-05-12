__precompile__()

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

    #sol = solveNetzwerk()
    #plotsol(sol)
end