Base.@kwdef mutable struct WPR_Param
    cv_H2O = 4182.0
    P0 = 1.0e5
    T0 = 293.15
    m_soll = 3.0
    Kp = 1.0e-4
end

Base.@kwdef mutable struct y_WPR
    Param::WPR_Param
    P::Number = Param.P0
    T::Number = Param.T0
end

Base.@kwdef mutable struct WPR_Knoten <: Wasser_Knoten
    #-- geänderte Parameter
    Param::WPR_Param

    #-- Zustandsvariablen
    y = y_WPR(Param=Param)
    
    #-- M-Matrix
    M::Array{Int} = [1; 1] 

    #-- zusätzeliche Infos
    Z::Dict
    
    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0   
    
    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end


function Knoten!(dy,k,knoten::WPR_Knoten,t)
    (; m_soll,Kp) = knoten.Param
    Masse = 0.1
    cv_H2O = 4182.0
    m = knoten.in[1].y.m
    dy[k] = Kp*(m_soll-m)
    dy[k+1] = knoten.sum_e/(Masse*cv_H2O) 
end