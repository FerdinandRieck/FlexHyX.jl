Base.@kwdef mutable struct GP0_Param
    P0 = 101325.0
    T0 = 293.15
end

Base.@kwdef mutable struct y_GP0
    Param::GP0_Param
    P::Number = Param.P0
    T::Number = Param.T0
end

Base.@kwdef mutable struct GP0_Knoten <: Gas_Knoten
    #-- geänderte Parameter
    Param::GP0_Param

    #-- Zustandsvariablen
    y = y_GP0(Param=Param)
    
    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Knoten!(dy,k,knoten::GP0_Knoten,t)
    #-- Parameter
    (; P0,T0) = knoten.Param
    #--

    P = knoten.y.P
    T = knoten.y.T
    
    dy[k] = P-P0
    dy[k+1] = T-T0
end