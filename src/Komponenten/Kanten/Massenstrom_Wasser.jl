Base.@kwdef mutable struct mWf2_Param
    m_dot = 1.0
    cv_H2O = 4182.0; #-- noch faglich
end

Base.@kwdef mutable struct y_mWf2
    m::Number = 0.0
    e::Number = 0.0
end

Base.@kwdef mutable struct mWf2_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWf2_Param

    #-- Zustandsvariablen
    y = y_mWf2()

    #-- Gasknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWf2_kante,t)
    #-- Parameter
    (; m_dot, cv_H2O) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    e = kante.y.e
    #--

    (; KL,KR,Z) = kante
    TL = KL.y.T
    TR = KR.y.T

    io = 1.0; if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end

    dy[k] = m - io*m_dot
    dy[k+1] = e - m*(TL-TR)
end