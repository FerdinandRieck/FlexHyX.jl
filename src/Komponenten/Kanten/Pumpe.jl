Base.@kwdef mutable struct mWPu_Param
    cv_H2O = 4182.0; #-- noch faglich
    m_max = 3.0 #-- max. Massenfluss = m_max
    a0 = 10.0 #-- max 10 bar Druckerhöhung
    Jac_init = true
end

Base.@kwdef mutable struct y_mWPu
    m::Number = 0.0
    e::Number = 0.0
end

Base.@kwdef mutable struct mWPu_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWPu_Param

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- Zustandsvariablen
    y = y_mWPu()

    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWPu_kante,t)
    #-- Parameter
    (; cv_H2O,a0,m_max,Jac_init) = kante.Param
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

    io = 1.0
    if Jac_init == true
        io = 1.0
    elseif (haskey(Z,"Schaltzeit_ein")==true) 
        io = ein(t,Z["Schaltzeit_ein"],Z["Schaltdauer"]) 
        io = io/length(Z["Schaltzeit_ein"])
    elseif (haskey(Z,"Schaltzeit_aus")==true) 
        io = aus(t,Z["Schaltzeit_aus"],Z["Schaltdauer"]) 
        io = io/length(Z["Schaltzeit_aus"])
    elseif (haskey(Z,"Schaltzeit")==true) 
        io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) 
    end
    a1 = a0 / m_max^2 
    dP = PR - PL; 
    dy[k] = m - io*sqrt(max(0.0,(a0-dP)/a1))
    dy[k+1] = e - 1.0e-6*cv_H2O*m*ifxaorb(m,TL,TR) 
end