Base.@kwdef mutable struct TMHR_Param
    T0 = 292.15
    Masse = 612.0
    c = 1006.7
end

Base.@kwdef mutable struct y_TMHR
    Param::TMHR_Param
    T::Number = Param.T0
end

Base.@kwdef mutable struct TMHR_Knoten <: Temp_Knoten
    #-- default Parameter
    Param::TMHR_Param

    #-- Zustandsvariable
    y = y_TMHR(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1] 

    #-- zusätzeliche Infos
    Z::Dict

    #-- Knotenbilanz
    sum_e::Number = 0.0  
    
    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::TMHR_Knoten,t)
    (; Masse, c) = knoten.Param

    eL = knoten.in[1].y.eL
    eR = knoten.in[1].y.eR

    e_in = eL-eR
    e_out = 0.0

    for i=1:length(knoten.out)
        e_out = e_out + knoten.out[i].y.e
    end
   
    dy[k] = (e_in - e_out)/(1e-6*Masse*c) 
end