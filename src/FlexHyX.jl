module FlexHyX
__precompile__()

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
    #=
    include("Funktionen/DGL.jl")
    include("Funktionen/Glättung.jl")
    include("Funktionen/Solve_Netzwerk.jl")
    include("Funktionen/Plot_Sol.jl")
    include("Funktionen/Read_Netz.jl")
    include("Funktionen/Inzidenz.jl")
    include("Funktionen/Datenstruktur.jl")
    include("Komponenten/Kanten/Batterie.jl")
    include("Komponenten/Kanten/Verbraucher.jl")
    include("Komponenten/Kanten/PV_Anlage.jl")
    include("Komponenten/Knoten/Kopplungsknoten_Strom.jl")
    include("Komponenten/Knoten/Potentialknoten.jl")
    =#
    #-- Funktionen exportieren
    export solveNetzwerk
    export plotsol

    #sol = solveNetzwerk()
    #plotsol(sol)
end

#Pkg.add(Pkg.PackageSpec(;name="DifferentialEquations", version="7.5.0"))