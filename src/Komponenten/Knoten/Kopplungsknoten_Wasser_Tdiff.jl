Base.@kwdef mutable struct WPTd_Param
    cv_H2O = 4182.0
    PW0 = 1.0e5
    T0 = 293.15
end

Base.@kwdef mutable struct y_WPTd
    Param::WPTd_Param
    P::Number = Param.PW0
    T::Number = Param.T0
end

Base.@kwdef mutable struct WPTd_Knoten <: Wasser_Knoten
    #-- geänderte Parameter
    Param::WPTd_Param

    #-- Zustandsvariablen
    y = y_WPTd(Param=Param)
    
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

function Knoten!(dy,k,knoten::WPTd_Knoten,t)
    (; cv_H2O) = knoten.Param

    Masse = 0.1
    dy[k] = sum_m(knoten.in,knoten.out)
    dy[k+1] = knoten.sum_e/(Masse*cv_H2O) 
end