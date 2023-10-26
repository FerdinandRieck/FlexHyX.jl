Base.@kwdef mutable struct WP0_Param
    PW0 = 1.013250
    T0 = 293.15
end

Base.@kwdef mutable struct y_WP0
    Param::WP0_Param
    P::Number = Param.PW0
    T::Number = Param.T0
end

Base.@kwdef mutable struct WP0_Knoten <: Wasser_Knoten
    #-- geänderte Parameter
    Param::WP0_Param

    #-- Zustandsvariablen
    y = y_WP0(Param=Param)
    
    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzliche Infos
    Z::Dict

    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::WP0_Knoten,t)
    #-- Parameter
    (; PW0,T0) = knoten.Param
    #--

    P = knoten.y.P
    T = knoten.y.T
    
    dy[k] = P-PW0
    dy[k+1] = T-T0
end