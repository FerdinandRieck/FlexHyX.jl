Base.@kwdef mutable struct mWPuE_Param
    cv_H2O = 4182.0; #-- noch faglich
end

Base.@kwdef mutable struct y_mWPuE
    m::Number = 0.0
end

Base.@kwdef mutable struct mWPuE_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWPuE_Param

    #-- Zustandsvariablen
    y = y_mWPuE()

    #-- Gasknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- M-Matrix
    M::Array{Int} = [0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWPuE_kante,t)
    #-- Parameter
    (; cv_H2O) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    #--

    (; KL,KR,Z) = kante
    PL = KL.y.P
    PR = KR.y.P
    TL = KL.y.T
    TR = KR.y.T

    io = 1.0; if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end
    a0 = 10; #-- max 10 bar Druckerhöhung
    m_max = 7; a1 = a0 / m_max^2 #-- max. Massenfluss = m_max
    dP = PR - PL; 
    dy[k] = m - io*sqrt(max(0.0,(a0-dP)/a1))
    e = 1e-6*cv_H2O*m*ifxaorb(m,TL,TR) 
    Z["e"] = e
end