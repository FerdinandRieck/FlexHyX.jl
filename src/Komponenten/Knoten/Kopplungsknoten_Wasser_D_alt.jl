Base.@kwdef mutable struct WPD_Param
    P0 = 1.0e5
end

Base.@kwdef mutable struct y_WPD
    Param::WPD_Param
    P::Number = Param.P0
end

Base.@kwdef mutable struct WPD_Knoten <: Wasser_Knoten
    #-- geänderte Parameter
    Param::WPD_Param

    #-- Zustandsvariablen
    y = y_WPD(Param=Param)
    
    #-- M-Matrix
    M::Array{Int} = [0] 

    #-- zusätzeliche Infos
    Z::Dict
    
    #-- Knotenbilanz
    sum_m::Number = 0.0
    
    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end


function Knoten!(dy,k,knoten::WPD_Knoten,t)
    dy[k] = knoten.sum_m
end