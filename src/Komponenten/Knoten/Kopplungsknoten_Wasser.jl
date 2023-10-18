Base.@kwdef mutable struct WP_Param
    cv_H2O = 4182.0
    PW0 = 1.0e5
    T0 = 293.15
end

Base.@kwdef mutable struct y_WP
    Param::WP_Param
    P::Number = Param.PW0
    T::Number = Param.T0
end

Base.@kwdef mutable struct WP_Knoten <: Wasser_Knoten
    #-- geänderte Parameter
    Param::WP_Param

    #-- Zustandsvariablen
    y = y_WP(Param=Param)
    
    #-- M-Matrix
    M::Array{Int} = [0; 1] 

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
    (; cv_H2O) = knoten.Param

    Masse = 0.1
    dy[k] = sum_m(knoten.in,knoten.out)
    dy[k+1] = knoten.sum_e/(1e-6*Masse*cv_H2O) 
end