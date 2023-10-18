Base.@kwdef mutable struct mWfPI_Param
    Kp = 0.001
    Ki = 0.0
    T_soll = 293.15
    ϕ0 = 0.0
    m_max = 1.0
    cv_H2O = 4182.0; #-- noch faglich
end

Base.@kwdef mutable struct y_mWfPI
    Param::mWfPI_Param
    m::Number = 0.0
    e::Number = 0.0
    T_err_dt::Number = 0.0
    ϕ::Number = Param.ϕ0
end

Base.@kwdef mutable struct mWfPI_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWfPI_Param

    #-- Zustandsvariablen
    y = y_mWfPI(Param=Param)

    #-- Wasserknoten links und rechts
    KL::Knoten
    KR::Knoten
    K = []

    #-- M-Matrix
    M::Array{Int} = [0; 0; 1; 1] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWfPI_kante,t)
    #-- Parameter
    (; T_soll, Kp, Ki, cv_H2O, m_max) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    e = kante.y.e
    T_err_dt = kante.y.T_err_dt
    ϕ = kante.y.ϕ 
    #--

    (; KL,KR,K,Z) = kante
    TL = KL.y.T
    TR = KR.y.T
    T = K.y.T

    io = 1.0; if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end
    ϕ_max=1; ϕ_min=0; ϕ = min(max(ϕ,ϕ_min),ϕ_max);

    dy[k] =  m - m_max*ϕ 
    dy[k+1] = e - 1e-6*cv_H2O*m*ifxaorb(m,TL,TR)
    dy[k+2] = T_soll-T
    dy[k+3] = Kp*(T_soll-T) + Ki*T_err_dt  
end