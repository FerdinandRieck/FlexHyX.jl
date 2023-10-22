Base.@kwdef mutable struct mWPuTd_Param
    cv_H2O = 4182.0; #-- noch faglich
    m_max = 3.0 #-- max. Massenfluss = m_max
    a0 = 10.0 #-- max 10 bar Druckerhöhung
end

Base.@kwdef mutable struct y_mWPuTd
    m::Number = 0.0
    e::Number = 0.0
end

Base.@kwdef mutable struct mWPuTd_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWPuTd_Param

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- Zustandsvariablen
    y = y_mWPuTd()

    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWPuTd_kante,t)
    #-- Parameter
    (; cv_H2O,m_max,a0) = kante.Param
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

    io = 1.0; if (haskey(Z,"Schaltzeit")==true) io = ein(t,Z["Schaltzeit"],Z["Schaltdauer"]) end
    io = io/length(Z["Schaltzeit"])
    a1 = a0 / m_max^2 
    dP = PR - PL; 

    dy[k] = m - io*sqrt(max(0.0,(a0-dP)/a1))
    dy[k+1] = e - cv_H2O*m*(TL-TR) 
end