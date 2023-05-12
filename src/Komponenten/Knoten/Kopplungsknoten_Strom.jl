Base.@kwdef mutable struct U_Param

end

Base.@kwdef mutable struct y_U
    U::Number = 0.0
end

Base.@kwdef mutable struct U_Knoten <: Strom_Knoten
    #-- geänderte Parameter
    Param::U_Param

    #-- Zustandsvariable
    y = y_U()

    #-- M-Matrix
    M::Array{Int} = [0] 

    #-- Jacobi Struktur
    J_fluss::Array{String} = ["sum_i"]

    #-- zusätzeliche Infos
    Z::Dict
end

function Knoten!(dy,k,sum_i,sum_m,sum_e,knoten::U_Knoten)
    dy[k] = sum_i
    nothing
end