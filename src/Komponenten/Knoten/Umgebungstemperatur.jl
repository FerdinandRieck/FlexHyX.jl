Base.@kwdef mutable struct T0_Param
    T0 = 293.15
end

Base.@kwdef mutable struct y_T0
    Param::T0_Param
    T::Number = Param.T0
end

Base.@kwdef mutable struct T0_Knoten <: Knoten
    #-- default Parameter
    Param::T0_Param

    #-- Zustandsvariable
    y = y_T0(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [0] 

    #-- zusätzeliche Infos
    Z::Dict

    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::T0_Knoten,t)
    (; T0) = knoten.Param
    T = knoten.y.T

    dy[k] = T-T0
end