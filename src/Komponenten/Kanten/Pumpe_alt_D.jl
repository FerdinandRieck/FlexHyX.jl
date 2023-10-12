Base.@kwdef mutable struct mWPuD_Param
    cv_H2O = 4182.0; #-- noch faglich
end

Base.@kwdef mutable struct y_mWPuD
    m::Number = 0.0
end

Base.@kwdef mutable struct mWPuD_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWPuD_Param

    #-- Zustandsvariablen
    y = y_mWPuD()

    #-- Gasknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- M-Matrix
    M::Array{Int} = [0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWPuD_kante,t)
    #-- Parameter
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    #--

    (; KL,KR,Z) = kante
    PL = KL.y.P
    PR = KR.y.P

    io = 1.0; if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end
    a0 = 1.0e5; #-- max 1 bar Druckerhöhung
    m_max = 1; a1 = a0 / m_max^2 #-- max. Massenfluss = m_max
    dP = PR - PL; 
    dy[k] = m - io*3.0 #sqrt(max(0.0,(a0-dP)/a1))
end