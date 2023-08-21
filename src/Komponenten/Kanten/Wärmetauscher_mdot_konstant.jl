Base.@kwdef mutable struct mWTaRK_Param
    nx = 1
    L = 0.286
    dx = L/max(nx,1/2)
    rho0 = 1000
    Dü = 0.02
    Da = 0.1
    A = pi/4*abs(Da^2-Dü^2)
    lamW = 0.6
    cv_H2O = 4182.0; #-- noch faglich
    Arho = A*rho0
    leit = lamW/(rho0*cv_H2O)
    kA = 4000.0
    WENO = true
    Richtung = "gleich"
    fluxTL = Array{Number}(undef, nx+1)
    fluxTR = Array{Number}(undef, nx+1)
    m_dot = 1.0
end

Base.@kwdef mutable struct y_mWTaRK
    m::Number = 0.0
    eL::Number = 0.0
    eR::Number = 0.0
    T
end


Base.@kwdef mutable struct mWTaRK_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWTaRK_Param

    #-- Knoten links und rechts
    KL::Knoten
    KR::Knoten

    #-- Zustandsvariablen
    y = y_mWTaRK(T = f(Vector(1:2:2*Param.nx), KL.y.T, KR.y.T, Param.nx))

    #-- Rohr aussen/innen
    R = 0

    #-- M-Matrix
    M::Array{Int} = [0; 0; 0; ones(Int,Param.nx)]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWTaRK_kante,t)
    #-- Parameter
    (; nx,dx,leit,Arho,Dü,cv_H2O,lamW,kA,WENO,Richtung,fluxTL,fluxTR) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    eL = kante.y.eL
    eR = kante.y.eR
    T = kante.y.T
    #--

    (; KL,KR,R) = kante
    TL = KL.y.T
    TR = KR.y.T
    mL = KL.in[1].y.m
    T_aussen = R.y.T

    if WENO == true
        recover_weno!(T,fluxTL,fluxTR) #??? hier vieleicht nur recover!(T), da bereits ubwind Diskertisierung mit ifxaorb(m[i],T[i]-fLT,fRT-T[i])??? 
    else
        recover!(T,fluxTL,fluxTR)
    end

    if (Richtung == "gegen") && nx > 1
        T_aussen = reverse(T_aussen)
    end

    #-- Rohr links
    dy[k] = m - mL #-- m
    dy[k+1] = eL - m*(TL-T[1])*cv_H2O  #- A/dx*2*lamW*(TL-T[1]) Wärmeübetragung weglassen
    #dy[k+1] = eL -(0.5*cv_H2O*(abs(m)*(TL-T[1])+m*(TL+T[1])) + A/dx*2*lamW*(TL-T[1])); #-- eL

    #-- Rohr rechts
    dy[k+2] = eR - m*(T[nx]-TR)*cv_H2O  #- A/dx*2*lamW*(T[nx]-TR) Wärmeübertragung weglassen
    #dy[k+2] = eR -(0.5*cv_H2O*(abs(m)*(T[end]-TR)+m*(T[end]+TR)) + A/dx*2*lamW*(T[end]-TR)) #-- eR

    #-- Rohr mitte
    fRT = TL #-- linke Randbedingung
    for i = 1:nx
        fLT = fRT
        if i == nx
            fRT = TR
        else 
            fRT = 0.5*(fluxTL[i+1]+fluxTR[i+1])
        end
        dy[k+i+2] = -1/Arho*m*ifxaorb(m,T[i]-fLT,fRT-T[i])*2/dx + leit*2/(dx^2)*(fLT-2*T[i]+fRT) - kA*pi*Dü/(Arho*cv_H2O)*(T[i]-T_aussen[i])  #-- T
    end
end