Base.@kwdef mutable struct WP_Param
    cv_H2O = 4182.0
    P0 = 1.0e5
    T0 = 293.15
end

Base.@kwdef mutable struct y_WP
    Param::WP_Param
    P::Number = Param.P0
    T::Number = Param.T0
end

Base.@kwdef mutable struct WP_Knoten <: Wasser_Knoten
    #-- geänderte Parameter
    Param::WP_Param

    #-- Zustandsvariablen
    y = y_WP(Param=Param)
    
    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzeliche Infos
    Z::Dict
    
    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0   
    
    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end


function Knoten!(dy,k,knoten::WP_Knoten,t)
    dy[k] = knoten.sum_m
    dy[k+1] = knoten.sum_e
end