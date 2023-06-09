Base.@kwdef mutable struct U0_Param
    U0 = 0.0
end

Base.@kwdef mutable struct y_U0
    Param::U0_Param
    U::Number = Param.U0
end

Base.@kwdef mutable struct U0_Knoten <: Strom_Knoten
    #-- default Parameter
    Param::U0_Param

    #-- Zustandsvariable (Kleinbuchstaben verwenden!)
    y = y_U0(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [0] 

    #-- Jacobi Struktur
    J::Int = 1

    #-- zusätzeliche Infos
    Z::Dict
end

function Knoten!(dy,k,knoten::U0_Knoten,t)
    (; U0) = knoten.Param
    U = knoten.y.U

    dy[k] = U-U0
end