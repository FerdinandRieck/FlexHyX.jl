Base.@kwdef mutable struct mWPu_Param
    cv_H2O = 4182.0; #-- noch faglich
end

Base.@kwdef mutable struct y_mWPu
    m::Number = 0.0
    e::Number = 0.0
end

Base.@kwdef mutable struct mWPu_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWPu_Param

    #-- Zustandsvariablen
    y = y_mWPu()

    #-- Gasknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWPu_kante,t)
    #-- Parameter
    (; cv_H2O) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    e = kante.y.e
    #--

    (; KL,KR,Z) = kante
    PL = KL.y.P
    PR = KR.y.P
    TL = KL.y.T
    TR = KR.y.T

    dy[k] = m - 15*(2 - 1.0e-5*(PR-PL))
    dy[k+1] = e - cv_H2O*m*ifxaorb(m,TL,TR)
end