Base.@kwdef mutable struct mWPu2_Param
    
end

Base.@kwdef mutable struct y_mWPu2
    m::Number = 0.0
    e::Number = 0.0
end

Base.@kwdef mutable struct mWPu2_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWPu2_Param

    #-- Zustandsvariablen
    y = y_mWPu2()

    #-- Gasknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWPu2_kante,t)
    #-- Parameter
    (; ) = kante.Param
    #--

    cv_H2O = 4182.0

    #-- Zustandsvariablen
    m = kante.y.m
    e = kante.y.e
    #--

    (; KL,KR,Z) = kante
    PL = KL.y.P
    PR = KR.y.P
    TL = KL.y.T
    TR = KR.y.T

    io = 1.0; if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end
    a0 = 1.0e5; #-- max 1 bar Druckerhöhung
    m_max = 1; a1 = a0 / m_max^2 #-- max. Massenfluss = m_max
    dP = PR - PL; 
    dy[k] = m - io*sqrt(max(0.0,(a0-dP)/a1))
    dy[k+1] = e - m*(TL-TR) 
end