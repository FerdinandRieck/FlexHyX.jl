module FlexHyX
    using DifferentialEquations, NLsolve, Plots
    using Symbolics
    using TerminalLoggers
    using Dates
    using LinearAlgebra
    using SparseArrays
    using SparseDiffTools
    using LinearSolve
    using Measures
    import JSON

    dir = dirname(@__FILE__)

    #-- Event-Funktion einfügen
    if ispath(pwd()*"/Events/")
        pfad = filter(contains(r".jl$"), readdir(pwd()*"/Events/";join=true))
        include.(pfad)
    end

    #-- Datenstruktur einfügen
    include("Datenstruktur.jl")

    #-- Knoten einfügen
    pfad = filter(contains(r".jl$"), readdir(dir*"/Komponenten/Knoten/";join=true))
    include.(pfad)

    #-- Kanten einfügen
    pfad = filter(contains(r".jl$"), readdir(dir*"/Komponenten/Kanten/";join=true))
    include.(pfad)
    
    #-- Funktionen einfügen
    pfad = filter(contains(r".jl$"), readdir(dir*"/Funktionen/";join=true))
    include.(pfad)

    #=
    #-- Typenhierarchie anzeigen
    using AbstractTrees
	AbstractTrees.children(x::Type) = subtypes(x)
    print_tree(flexhyx)
    =#

    #-- Funktionen exportieren
    export solveNetz, plotSol, plotNetz
 end