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

    #-- zusätzeliche Infos
    Z::Dict

    #-- Knotenbilanz
    sum_i::Number = 0.0  

    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::U_Knoten,t)
    dy[k] = knoten.sum_i
    nothing
end

