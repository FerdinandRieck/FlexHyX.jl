Base.@kwdef mutable struct mWVentil_Param
    KV = 7.0*1.0e-8
    P_soll = 1.2
    alpha = 2.0*1.0e7
    cv_H2O = 4182.0; #-- noch faglich
    A0 = 0.5
end

Base.@kwdef mutable struct y_mWVentil
    Param::mWVentil_Param
    m::Number = 0.0
    e::Number = 0.0
    A::Number = Param.A0
end

Base.@kwdef mutable struct mWVentil_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWVentil_Param

    #-- Zustandsvariablen
    y = y_mWVentil(Param=Param)

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- M-Matrix
    M::Array{Int} = [0; 0; 1] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWVentil_kante,t)
    #-- Parameter
    (; KV,P_soll,alpha,cv_H2O) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    e = kante.y.e
    A = kante.y.A;
    #--

    (; KL,KR,Z) = kante
    PL = KL.y.P
    PR = KR.y.P
    TL = KL.y.T
    TR = KR.y.T

    io = 1.0; if get(Z,"Schaltzeit",0)!=0 io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end
    P = PL-PR;  
    A_max=1; A_min=0; A = min(max(A,A_min),A_max);
    diff = P_soll*1.0e5 - PR;

    dy[k] = m - KV*A*P/(sqrt(abs(P)+100));
    dy[k+1] = e - cv_H2O*m*ifxaorb(m,TL,TR)
    dy[k+2] = diff/alpha*ifxaorb(diff,A_max-A,A-A_min) - (1-io)*A
end