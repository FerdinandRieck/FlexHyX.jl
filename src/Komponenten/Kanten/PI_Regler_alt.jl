Base.@kwdef mutable struct mPI_Param
    Kp = 0.001
    Ki = 0.005
    T_soll = 293.15
    cv_H2O = 4182.0; #-- noch faglich
end

Base.@kwdef mutable struct y_mPI
    m::Number = 0.0
    T_err_dt::Number = 0.0
    e::Number = 0.0
end

Base.@kwdef mutable struct mPI_kante <: Wasser_Kante
    #-- default Parameter
    Param::mPI_Param

    #-- Zustandsvariablen
    y = y_mPI()

    #-- Gasknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten
    K = []

    #-- M-Matrix
    M::Array{Int} = [1; 1; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mPI_kante,t)
    #-- Parameter
    (; T_soll, Kp, Ki, cv_H2O) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    T_err_dt = kante.y.T_err_dt
    e = kante.y.e
    #--

    (; KL,KR,K,Z) = kante
    TL = KL.y.T
    TR = KR.y.T
    T = K.y.T

    io = 1.0; if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end

    u = Kp*(T_soll-T) + Ki*T_err_dt 

    dy[k] = u 
    dy[k+1] = T_soll-T
    dy[k+2] = e - 1e-6*cv_H2O*m*ifxaorb(m,TL,TR)
end