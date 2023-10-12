Base.@kwdef mutable struct mWPuR_Param
    cv_H2O = 4182.0; #-- noch faglich
    m = 0.0
    m_max = 3.0
    P_Soll = 1.51
    ϕ0 = 0.0
    Kp = 0.01
end

Base.@kwdef mutable struct y_mWPuR
    Param::mWPuR_Param
    m::Number = Param.m
    e::Number = 0.0
    ϕ::Number = Param.ϕ0
end

Base.@kwdef mutable struct mWPuR_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWPuR_Param

    #-- Zustandsvariablen
    y = y_mWPuR(Param=Param)

    #-- Gasknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- M-Matrix
    M::Array{Int} = [0; 0; 1] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWPuR_kante,t)
    #-- Parameter
    (; cv_H2O,P_Soll,Kp,m_max) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    e = kante.y.e
    ϕ = kante.y.ϕ
    #--

    (; KL,KR,Z) = kante
    PL = KL.y.P
    PR = KR.y.P
    TL = KL.y.T
    TR = KR.y.T

    io = 1.0; if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end
    P_Soll = P_Soll*1.0e5
    ϕ_max = 1.0; ϕ_min = 0.0; 
    ϕ = min(max(ϕ,ϕ_min),ϕ_max)

    dy[k] = m - ϕ*m_max # io*sqrt(max(0.0,(a0-dP)/a1))
    dy[k+1] = e - cv_H2O*m*ifxaorb(m,TL,TR)
    dy[k+2] = Kp*(P_Soll-PR)*ifxaorb((P_Soll-PR),ϕ_max-ϕ,ϕ-ϕ_min) - (1-io)*ϕ 
end