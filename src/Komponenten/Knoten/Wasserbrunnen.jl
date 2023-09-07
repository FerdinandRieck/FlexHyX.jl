Base.@kwdef struct WPB2_Param
    PA = 101325
    P0 = 1.0e5
    T0 = 293.15
    D = 1.0
    A = pi*(D/2)^2
    cv_H2O = 4182; #-- noch faglich
end

Base.@kwdef mutable struct y_WPB2
    Param::WPB2_Param
    M::Number = Param.P0*Param.A/9.81
    T::Number = Param.T0
    P::Number = Param.P0
    P_boden::Number = Param.P0
end

Base.@kwdef mutable struct WPB2_Knoten <: Wasser_Knoten
    #-- default Parameter
    Param::WPB2_Param

    #-- Zustandsvariablen 
    y = y_WPB2(Param=Param)

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

function Knoten!(dy,k,knoten::WPB2_Knoten,t)
    #-- Parameter
    (; A,cv_H2O,PA) = knoten.Param
    #--

    (; M, P, P_boden) = knoten.y

    dy[k] = knoten.sum_m
    dy[k+1] = knoten.sum_e/(M*cv_H2O)
    dy[k+2] = P-PA
    dy[k+3] = P_boden-M*9.81/A
end