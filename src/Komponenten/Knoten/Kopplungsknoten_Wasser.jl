Base.@kwdef mutable struct WP2_Param

end

Base.@kwdef mutable struct y_WP2
    P::Number = 0.0
    T::Number = 293.15
end

Base.@kwdef mutable struct WP2_Knoten <: Wasser_Knoten
    #-- geänderte Parameter
    Param::WP2_Param

    #-- Zustandsvariablen
    y = y_WP2()
    
    #-- M-Matrix
    M::Array{Int} = [0; 1] 

    #-- zusätzeliche Infos
    Z::Dict
    
    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0   
    
    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end


function Knoten!(dy,k,knoten::WP2_Knoten,t)
    Masse = 0.1
    dy[k] = sum_m(knoten.in,knoten.out)
    dy[k+1] =  knoten.sum_e/Masse 
end