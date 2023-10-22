Base.@kwdef mutable struct GP_Param

end

Base.@kwdef mutable struct y_GP
    P::Number = 0.0
    T::Number = 293.15
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
    
    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0   
    
    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end


function Knoten!(dy,k,knoten::GP_Knoten,t)
    dy[k] = knoten.sum_m
    dy[k+1] = knoten.sum_e
end