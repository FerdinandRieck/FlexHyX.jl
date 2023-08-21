Base.@kwdef mutable struct mPID_Param
    Kp = 0.001
    Ki = 0.005
    Kd = 0.001
    T_soll = 293.15
    cv_H2O = 4182.0; #-- noch faglich
end

Base.@kwdef mutable struct y_mPID
    Edot::Number = 0.0
    m::Number = 0.0
    Edt::Number = 3.0
    e::Number = 0.0
end

Base.@kwdef mutable struct mPID_kante <: Wasser_Kante
    #-- default Parameter
    Param::mPID_Param

    #-- Zustandsvariablen
    y = y_mPID()

    #-- Gasknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten
    K = []

    #-- M-Matrix
    M::Array{Int} = [1; 1; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mPID_kante,t)
    #-- Parameter
    (; T_soll, Kp, Ki, Kd, cv_H2O) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    T_err_dt = kante.y.T_err_dt
    e = kante.y.e
    #--

    (; KL,KR,K) = kante
    TL = KL.y.T
    TR = KR.y.T
    T = K.y.T

    E = T - T_soll
    dy[k] = m
    dy[k+1] = - Kp*E - Ki*Edt - Kd*m
    dy[k+2] = E
    dy[k+2] = e - m*(TL-TR)*cv_H2O 
end


