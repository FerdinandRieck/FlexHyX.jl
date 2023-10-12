Base.@kwdef mutable struct mWPuP_Param
    cv_H2O = 4182.0; #-- noch faglich
    m = 0.0
    m_max = 3.0
    P_Soll = 1.51
    Kp = 0.01
end

Base.@kwdef mutable struct y_mWPuP
    Param::mWPuP_Param
    m::Number = Param.m
    e::Number = 0.0
end

Base.@kwdef mutable struct mWPuP_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWPuP_Param

    #-- Zustandsvariablen
    y = y_mWPuP(Param=Param)

    #-- Gasknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- M-Matrix
    M::Array{Int} = [1; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWPuP_kante,t)
    #-- Parameter
    (; cv_H2O,P_Soll,Kp) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    e = kante.y.e
    #--

    (; KL,KR) = kante
    PL = KL.y.P
    PR = KR.y.P
    TL = KL.y.T
    TR = KR.y.T

    P_Soll = P_Soll*1.0e5

    dy[k] = Kp*(P_Soll-PR)
    dy[k+1] = e - cv_H2O*m*ifxaorb(m,TL,TR)
end