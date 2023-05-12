Base.@kwdef mutable struct TM_Param
    T0 = 293.15
    Masse = 10
    c = 896
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

    #-- Jacobi Struktur

    #-- zusätzeliche Infos
    Z::Dict
end

function Knoten!(dy,k,sum_i,sum_m,sum_e,knoten::TM_Knoten)
    (; Masse, c) = knoten.Param

    dy[k] = sum_e/(c * Masse)
end