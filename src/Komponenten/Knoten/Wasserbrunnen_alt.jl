Base.@kwdef struct WPB_Param
    PA = 101325
    P0 = 1.0e5
    T0 = 293.15
    D = 1.0
    A = pi*(D/2)^2
    cv_H2O = 4182; #-- noch faglich
end

Base.@kwdef mutable struct y_WPB
    Param::WPB_Param
    M::Number = Param.P0*Param.A/9.81
    MT::Number = M*Param.T0 
    P::Number = Param.PA
    T::Number = Param.T0
    P_boden::Number = Param.P0
end

Base.@kwdef mutable struct WPB_Knoten <: Wasser_Knoten
    #-- default Parameter
    Param::WPB_Param

    #-- Zustandsvariablen 
    y = y_WPB(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1; 1; 0; 0; 0] 

    #-- zusätzeliche Infos
    Z::Dict

    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0
end

function Knoten!(dy,k,knoten::WPB_Knoten,t)
    #-- Parameter
    (; A,cv_H2O,PA) = knoten.Param
    #--

    (; M, MT, P, T, P_boden) = knoten.y

    dy[k] = knoten.sum_m
    dy[k+1] = knoten.sum_e/cv_H2O
    dy[k+2] = P-PA
    dy[k+3] = T-MT/M
    dy[k+4] = P_boden-M*9.81/A
end