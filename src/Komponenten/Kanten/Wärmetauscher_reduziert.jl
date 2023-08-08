Base.@kwdef mutable struct mWTaR2_Param
    nx = 1
    L = 1.0
    dx = L/max(nx,1/2)
    dx2 = L/max(1,1/2)
    D = 0.1
    A = pi*(D/2)^2
    g = 9.81
    a = 1414
    a2 = a^2
    mu = 1.0e-3
    rho0 = 1000.0
    lamW = 0.6
    cv_H2O = 4182.0; #-- noch faglich
    Arho = A*rho0
    leit = lamW/(rho0*cv_H2O)
    K = 1e-5 #-- Rauheit
    phi = 0.0 #-- Neigungswinkel
    kA = 380.0
    WENO = true
    Richtung = "gleich"
    fluxTL = Array{Number}(undef, nx+1)
    fluxTR = Array{Number}(undef, nx+1)
end

Base.@kwdef mutable struct y_mWTaR2
    mL::Number = 0.0
    eL::Number = 0.0
    P
    _m::Number = 0.0
    T
    mR::Number = 0.0
    eR::Number = 0.0
end


Base.@kwdef mutable struct mWTaR2_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWTaR2_Param

    #-- Knoten links und rechts
    KL::Knoten
    KR::Knoten

    #-- Zustandsvariablen
    y = y_mWTaR2(P = KL.y.P + 0.5*(KR.y.P-KL.y.P), 
                 T = f(Vector(1:2:2*Param.nx), KL.y.T, KR.y.T, Param.nx)
                )

    #-- Rohr aussen
    RA = 0

    #-- M-Matrix
    M::Array{Int} = [1; 0; 1; 1; ones(Int,Param.nx); 1; 0]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWTaR2_kante,t)
    #-- Parameter
    (; nx,dx,dx2,a2,leit,Arho,rho0,A,D,cv_H2O,mu,K,lamW,phi,g,kA,WENO,Richtung,fluxTL,fluxTR) = kante.Param
    #--

    #-- Zustandsvariablen
    mL = kante.y.mL
    eL = kante.y.eL
    P = kante.y.P
    m = kante.y._m
    T = kante.y.T
    mR = kante.y.mR
    eR = kante.y.eR
    #--

    (; KL,KR,RA) = kante
    PL = KL.y.P
    TL = KL.y.T
    PR = KR.y.P
    TR = KR.y.T
    T_aussen = RA.y.T

    if WENO == true
        recover_weno!(T,fluxTL,fluxTR) #??? hier vieleicht nur recover!(T), da bereits ubwind Diskertisierung mit ifxaorb(m[i],T[i]-fLT,fRT-T[i])??? 
    else
        recover!(T,fluxTL,fluxTR)
    end

    if (Richtung == "gegen") && nx > 1
        T_aussen = reverse(T_aussen)
    end

    #-- Rohr links
    dy[k] = -(m^2-mL^2)*2/(dx2*Arho) - A*(P-PL)*2/dx2 - lambda(mL,D,A,mu,K)/(2*D*Arho)*abs(mL)*mL - g*Arho*sin(phi)  #-- mL
    dy[k+1] = eL -(0.5*cv_H2O*(abs(mL)*(TL-T[1])+mL*(TL+T[1])) + A/dx2*2*lamW*(TL-T[1]))  #-- eL
    dy[k+2] = -a2/A*(mR-mL)/dx2  #-- P 
    dy[k+3] = -(mR^2-mL^2)/(dx2*Arho) - A*(PR-PL)/dx2 - lambda(m,D,A,mu,K)/(2*D*Arho)*abs(m)*m - g*Arho*sin(phi)  #-- m #!!! Nicht klar ob Winkel Angegeben werden darf, wegen GL.88 (Herleitung aus Text Rohr (Wasser))
    #-- Rohr mitte
    fRT = TL #-- linke Randbedingung
    for i = 1:nx
        fLT = fRT
        if i == nx
            fRT = TR
        else 
            fRT = 0.5*(fluxTL[i+1]+fluxTR[i+1])
        end
        dy[k+i+3] = -1/Arho*m*ifxaorb(m,T[i]-fLT,fRT-T[i])*2/dx + leit*2/(dx^2)*(fLT-2*T[i]+fRT) - 4*kA/(rho0*D*cv_H2O)*(T[i]-T_aussen[i])  #-- T
    end

    #-- Rohr rechts
    dy[k+nx+4] = -(mR^2-m^2)*2/(dx2*Arho) - A*(PR-P)*2/dx2 - lambda(mR,D,A,mu,K)/(2*D*Arho)*abs(mR)*mR - g*Arho*sin(phi)  #-- mR
    dy[k+nx+5] = eR -(0.5*cv_H2O*(abs(mR)*(T[end]-TR)+mR*(T[end]+TR)) + A/dx2*2*lamW*(T[end]-TR)) #-- eR
end