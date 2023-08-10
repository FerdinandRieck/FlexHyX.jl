Base.@kwdef mutable struct mWf_Param
    m_dot = 1.0
end

Base.@kwdef mutable struct y_mWf
    m::Number = 0.0
    e::Number = 0.0
end

Base.@kwdef mutable struct mWf_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWf_Param

    #-- Zustandsvariablen
    y = y_mWf()

    #-- Gasknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWf_kante,t)
    #-- Parameter
    (; m_dot) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    e = kante.y.e
    #--

    (; KL,KR,Z) = kante
    TL = KL.y.T
    TR = KR.y.T

    io = 1.0; if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end

    dy[k] = m - io*m_dot #- (1 - io)*0.01
    dy[k+1] = e - m*(TL-TR) 
end