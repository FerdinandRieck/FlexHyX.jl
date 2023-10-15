Base.@kwdef mutable struct T_Param
    cv_H2O = 4182.0
    T0 = 293.15
end

Base.@kwdef mutable struct y_T
    Param::T_Param
    T::Number = Param.T0
end

Base.@kwdef mutable struct T_Knoten <: Wasser_Knoten
    #-- geänderte Parameter
    Param::T_Param

    #-- Zustandsvariablen
    y = y_T(Param=Param)
    
    #-- M-Matrix
    M::Array{Int} = [1] 

    #-- zusätzeliche Infos
    Z::Dict
    
    #-- Knotenbilanz
    sum_e::Number = 0.0    

    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end


function Knoten!(dy,k,knoten::T_Knoten,t)
    (; cv_H2O) = knoten.Param
    Masse = 0.1
    dy[k] =  knoten.sum_e/(1e-6*Masse*cv_H2O)
end