module FlexHyX
    using DifferentialEquations, NLsolve, Plots
    using TerminalLoggers
    using Dates
    using LinearAlgebra
    using SparseArrays
    import Pkg
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

    #sol = solveNetzwerk(dir)
    #plotsol(sol,[38,46],28) # i_out -> iDiode 
    #plotsol(sol,[37,39],28) # iBZ   -> i_in
    #plotsol(sol,[40,35,42],28) # Summe m Knoten 14 
 end

#=
Aktuelle Bugs:
- Bei Änderung in Eventsfunktion wird erst beim zweiten Durchlauf Änderung übernommen
- Problem mit Events und normelne Run button
- leichte unterschiede der Eventzeitpunkte von Matlab und JULIA mit Zeitreihe
=#