Base.@kwdef mutable struct mWfP_Param
    Kp = 0.01
    cv_H2O = 4182.0; #-- noch faglich
    ϕ0 = 0.0
    T_soll = 293.15
    m_max = 1.0
end

Base.@kwdef mutable struct y_mWfP
    Param::mWfP_Param
    m::Number = 0.0
    e::Number = 0.0
    ϕ::Number = Param.ϕ0
end

Base.@kwdef mutable struct mWfP_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWfP_Param

    #-- Zustandsvariablen
    y = y_mWfP(Param=Param)

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten
    K = []

    #-- M-Matrix
    M::Array{Int} = [0; 0; 1] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWfP_kante,t)
    #-- Parameter
    (; T_soll,Kp,cv_H2O,m_max) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    e = kante.y.e
    ϕ = kante.y.ϕ;
    #--

    (; KL,KR,K,Z) = kante
    TL = KL.y.T
    TR = KR.y.T
    T = K.y.T

    #T_soll = T_soll + 273.15

    io = 1.0; if get(Z,"Schaltzeit",0)!=0 io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end 
    ϕ_max=1; ϕ_min=0; ϕ = min(max(ϕ,ϕ_min),ϕ_max);

    dy[k] = m - m_max*ϕ
    dy[k+1] = e - 1e-6*cv_H2O*m*ifxaorb(m,TL,TR)
    dy[k+2] = Kp*(T_soll-T)
end
