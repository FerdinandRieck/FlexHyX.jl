Base.@kwdef struct WPSPD_Param
    P0 = 1.0e5
    D = 1.0
    A = pi*(D/2)^2
    cv_H2O = 4182; #-- noch faglich
end

Base.@kwdef mutable struct y_WPSPD
    Param::WPSPD_Param
    M::Number = Param.P0*Param.A/9.81
    P::Number = Param.P0
end

Base.@kwdef mutable struct WPSPD_Knoten <: Wasser_Knoten
    #-- default Parameter
    Param::WPSPD_Param

    #-- Zustandsvariablen 
    y = y_WPSPD(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1; 0] 

    #-- zusätzeliche Infos
    Z::Dict

    #-- Knotenbilanz
    sum_m::Number = 0.0

    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::WPSPD_Knoten,t)
    #-- Parameter
    (; A,cv_H2O) = knoten.Param
    #--

    (; M, P) = knoten.y

    dy[k] = knoten.sum_m
    dy[k+1] = P-M*9.81/A
end
