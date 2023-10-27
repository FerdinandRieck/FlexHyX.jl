Base.@kwdef mutable struct mWf_Param
    m_dot = 1.0
    cv_H2O = 4182.0; #-- noch faglich
    Jac_init = true
end

Base.@kwdef mutable struct y_mWf
    m::Number = 0.0
    e::Number = 0.0
end

Base.@kwdef mutable struct mWf_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWf_Param

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- Zustandsvariablen
    y = y_mWf()

    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWf_kante,t)
    #-- Parameter
    (; m_dot,cv_H2O,Jac_init) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    e = kante.y.e
    #--

    (; KL,KR,Z) = kante
    TL = KL.y.T
    TR = KR.y.T

    io = 1.0
    if Jac_init == true
        io = 1.0
    elseif (haskey(Z,"Schaltzeit")==true) 
        io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) 
    end

    dy[k] = m - io*m_dot
    dy[k+1] = e - 1e-6*cv_H2O*m*ifxaorb(m,TL,TR)
end