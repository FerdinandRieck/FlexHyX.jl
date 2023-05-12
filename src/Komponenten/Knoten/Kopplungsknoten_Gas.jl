Base.@kwdef mutable struct GP_Param

end

Base.@kwdef mutable struct y_GP
    P::Number = 0.0
    T::Number = 0.0
end

Base.@kwdef mutable struct GP_Knoten <: Gas_Knoten
    #-- geänderte Parameter
    Param::GP_Param

    #-- Zustandsvariablen
    y = y_GP()
    
    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzeliche Infos
    Z::Dict
end


function Knoten!(dy,k,sum_i,sum_m,sum_e,knoten::GP_Knoten)
    dy[k] = sum_m
    dy[k+1] = sum_e
end