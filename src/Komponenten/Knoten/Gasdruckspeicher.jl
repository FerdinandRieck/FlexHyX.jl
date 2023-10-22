Base.@kwdef mutable struct GPSP_Param
    P0 = 144210.0
    T0 = 293.15
    V = 1.0
    Rs = 4124.2
    cv = 1.01798*1.0e4
end

Base.@kwdef mutable struct y_GPSP
    Param::GPSP_Param
    M::Number = Param.P0*Param.V/(Param.T0*Param.Rs)
    MT::Number = M*Param.T0 
    P::Number = Param.P0
    T::Number = Param.T0
end

Base.@kwdef mutable struct GPSP_Knoten <: Gas_Knoten
    #-- default Parameter
    Param::GPSP_Param

    #-- Zustandsvariablen 
    y = y_GPSP(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1; 1; 0; 0] 

    #-- zusätzeliche Infos
    Z::Dict

    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0

    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::GPSP_Knoten,t)
    #-- Parameter
    (; V,Rs,cv) = knoten.Param
    #--
    (; M, MT, P, T) = knoten.y
    dy[k] = knoten.sum_m
    dy[k+1] = knoten.sum_e/cv
    dy[k+2] = P-Rs*MT/V
    dy[k+3] = T-MT/max(M,1.0e-9)
end