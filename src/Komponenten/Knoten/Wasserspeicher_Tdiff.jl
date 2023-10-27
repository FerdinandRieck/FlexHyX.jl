Base.@kwdef mutable struct WPSPTd_Param
    PW0 = 1.0
    T0 = 293.15
    D = 1.0
    A = pi*(D/2)^2
    cv_H2O = 4182; #-- noch faglich
end

Base.@kwdef mutable struct y_WPSPTd
    Param::WPSPTd_Param
    M::Number = Param.PW0*1e5*Param.A/9.81
    T::Number = Param.T0
    P::Number = Param.PW0
end

Base.@kwdef mutable struct WPSPTd_Knoten <: Wasser_Knoten
    #-- default Parameter
    Param::WPSPTd_Param

    #-- Zustandsvariablen 
    y = y_WPSPTd(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1; 1; 0] 

    #-- zusätzeliche Infos
    Z::Dict

    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0

    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::WPSPTd_Knoten,t)
    #-- Parameter
    (; A,cv_H2O) = knoten.Param
    #--

    (; M, P) = knoten.y

    dy[k] =  knoten.sum_m
    dy[k+1] = knoten.sum_e/(M*cv_H2O)
    dy[k+2] = P-M*9.81/A*1e-5
end
