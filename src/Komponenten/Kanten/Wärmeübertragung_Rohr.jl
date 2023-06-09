Base.@kwdef mutable struct mWTR_Param
    L = 1.0
    dx = L/max(1,1/2)
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
    kA = 25.0
    T_aussen = 293.15
end

Base.@kwdef mutable struct y_mWTR
    mL::Number = 0.0
    eL::Number = 0.0
    P::Number = 1.0 
    _m::Number = 0.0
    T::Number = 3.0
    mR::Number = 0.0
    eR::Number = 0.0
end

Base.@kwdef mutable struct mWTR_kante <: Temp_Kante
    #-- default Parameter
    Param::mWTR_Param

    #-- Zustandsvariablen
    y = y_mWTR()

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten  
    KR::Wasser_Knoten  

    #-- M-Matrix
    M::Array{Int} = [1; 0; 1; 1; 1; 1; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWTR_kante,t)
    #-- Parameter
    (; a2,A,D,dx,Arho,leit,K,mu,phi,g,cv_H2O,lamW,kA,T_aussen) = kante.Param
    #--

    #-- Zustandsvariable
    mL = kante.y.mL
    eL = kante.y.eL
    mR = kante.y.mR
    eR = kante.y.eR
    P = kante.y.P
    m = kante.y._m
    T = kante.y.T
    #--

    (; KL,KR) = kante
    PL = KL.y.P
    TL = KL.y.T
    PR = KR.y.P
    TR = KR.y.T

    dy[k] = -(m^2-mL^2)*2/(dx*Arho) - A*(P-PL)*2/dx - lambda(mL,D,A,mu,K)/(2*D*Arho)*abs(mL)*mL - g*Arho*sin(phi) #-- mL
    dy[k+1] = eL -(0.5*cv_H2O*(abs(mL)*(TL-T)+mL*(TL+T)) + A/dx*2*lamW*(TL-T)); #-- eL
    dy[k+2] = -a2/A*(mR-mL)/dx  #-- P
    dy[k+3] = -(mR^2-mL^2)/(dx*Arho) - A*(PR-PL)/dx - lambda(m,D,A,mu,K)/(2*D*Arho)*abs(m)*m - g*Arho*sin(phi)  #-- m
    dy[k+4] = -1/Arho*m*ifxaorb(m,T-TL,TR-T)*2/dx + leit*2/(dx^2)*(TL-2*T+TR) - kA/(Arho*D)*(T-T_aussen)  #-- T
    dy[k+5] = -(mR^2-m^2)*2/(dx*Arho) - A*(PR-P)*2/dx - lambda(mR,D,A,mu,K)/(2*D*Arho)*abs(mR)*mR - g*Arho*sin(phi) #-- mR
    dy[k+6] = eR -(0.5*cv_H2O*(abs(mR)*(T-TR)+mR*(T+TR)) + A/dx*2*lamW*(T-TR)) #-- eR
end