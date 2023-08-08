Base.@kwdef mutable struct WP_Param

end

Base.@kwdef mutable struct y_WP
    P::Number = 0.0
    T::Number = 293.15
end

Base.@kwdef mutable struct WP_Knoten <: Wasser_Knoten
    #-- geänderte Parameter
    Param::WP_Param

    #-- Zustandsvariablen
    y = y_WP()
    
    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzeliche Infos
    Z::Dict
    
    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0    
end


function Knoten!(dy,k,knoten::WP_Knoten,t)
    dy[k] = knoten.sum_m
    dy[k+1] = knoten.sum_e
end