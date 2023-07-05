module FlexHyX
    using DifferentialEquations, NLsolve, Plots
    using Symbolics
    using TerminalLoggers
    using Dates
    using LinearAlgebra
    using SparseArrays
    import JSON

    dir = dirname(@__FILE__)

    #-- Event-Funktion einfügen
    if ispath(pwd()*"/Events/")
        pfad = filter(contains(r".jl$"), readdir(pwd()*"/Events/";join=true))
        include.(pfad)
    end
    
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
    print_tree(flexhyx)
    =#

    #-- Funktionen exportieren
    export solveNetzwerk
    export plotSol
    export plotNetz
 end

#=
Aktuelle Bugs:
- leichte unterschiede der Eventzeitpunkte von Matlab und JULIA mit Zeitreihe
=#

#=
große matrixoperation / sum_i Speicher nicht allokiert:
BenchmarkTools.Trial: 1 sample with 1 evaluation.
 Single result which took 9.941 s (5.38% GC) to evaluate,
 with a memory estimate of 2.62 GiB, over 94404818 allocations.

große matrixoperation / sum_i Speicher allokiert:
BenchmarkTools.Trial: 1 sample with 1 evaluation.
 Single result which took 10.184 s (5.54% GC) to evaluate,
 with a memory estimate of 2.62 GiB, over 94409015 allocations.

 kleine Matrixoperation / sum_i Speicher allokiert:
 BenchmarkTools.Trial: 1 sample with 1 evaluation.
 Single result which took 5.796 s (2.12% GC) to evaluate,
 with a memory estimate of 671.57 MiB, over 10987182 allocations.
=#