Base.@kwdef mutable struct TM_Param
    T0 = 293.15
    Masse = 10.0
    c = 896.0
end

Base.@kwdef mutable struct y_TM
    Param::TM_Param
    T::Number = Param.T0
end

Base.@kwdef mutable struct TM_Knoten <: Temp_Knoten
    #-- default Parameter
    Param::TM_Param

    #-- Zustandsvariable
    y = y_TM(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1] 

    #-- zusätzeliche Infos
    Z::Dict

    #-- Knotenbilanz
    sum_e::Number = 0.0    
end

function Knoten!(dy,k,knoten::TM_Knoten,t)
    (; Masse, c) = knoten.Param

    dy[k] = knoten.sum_e/(c )
end