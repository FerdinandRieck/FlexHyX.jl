Base.@kwdef mutable struct mWRo_Param
    nx = 1.0
    L = 1.0
    dx = L/max(nx,1/2)
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
end

#-- Rohrränder ---------------------------------
Base.@kwdef mutable struct y_mWRo_randL
    mL::Number = 0.0
    eL::Number = 0.0
end

Base.@kwdef mutable struct mWRo_kante_randL <: Wasser_Kante
    #-- default Parameter
    Param::mWRo_Param

    #-- Zustandsvariablen
    y = y_mWRo_randL()

    #-- Wasserknoten links
    KL::Wasser_Knoten

    #-- Rohrabschnitt links 
    RL = 0 

    #-- M-Matrix
    M::Array{Int} = [1; 0]

    #-- zusätzliche Infos
    Z::Dict
end

Base.@kwdef mutable struct y_mWRo_randR
    mR::Number = 0.0
    eR::Number = 0.0
end

Base.@kwdef mutable struct mWRo_kante_randR <: Wasser_Kante
    #-- default Parameter
    Param::mWRo_Param

    #-- Zustandsvariablen
    y = y_mWRo_randR()

    #-- Wasserknoten rechts
    KR::Wasser_Knoten

    #-- Rohrabschnitt rechts 
    RR = 0

    #-- M-Matrix
    M::Array{Int} = [1; 0]

    #-- zusätzliche Infos
    Z::Dict
end
#-----------------------------------------------

#-- Rohrmitte ----------------------------------
Base.@kwdef mutable struct y_mWRo_mitte
    P::Number = 1.0 #???Warum diese Anfangswerte in MATLAB???
    m::Number = 0.0
    T::Number = 3.0
end

Base.@kwdef mutable struct mWRo_kante_mitte <: Wasser_Kante
    #-- default Parameter
    Param::mWRo_Param

    #-- Zustandsvariablen
    y = y_mWRo_mitte()

    #-- Rohrabschnitt links und rechts
    RL = 0  
    RR = 0 

    #-- M-Matrix
    M::Array{Int} = [1; 1; 1]

    #-- zusätzliche Infos
    Z::Dict
end
#-----------------------------------------------

function Kante!(dy,k,kante::mWRo_kante_randL,t)
    #-- Parameter
    (; dx,Arho,A,D,cv_H2O,mu,K,lamW,phi,g) = kante.Param
    #--

    #-- Zustandsvariablen
    mL = kante.y.mL
    eL = kante.y.eL
    #--

    (; KL,RL) = kante
    PL = KL.y.P
    TL = KL.y.T
    m1 = RL.y.m
    P1 = RL.y.P
    T1 = RL.y.T

    lam = lamda(mL,D,A,mu,K)

    dy[k] = -(m1^2-mL^2)*2/(dx*Arho) - A*(P1-PL)*2/dx - lam/(2*D*Arho)*abs(mL)*mL - g*Arho*sin(phi); #-- mL
    dy[k+1] = eL -(0.5*cv_H2O*(abs(mL)*(TL-T1)+mL*(TL+T1)) + A/dx*2*lamW*(TL-T1)); #-- eL
end

function Kante!(dy,k,kante::mWRo_kante_randR,t)
    #-- Parameter
    (; dx,Arho,A,D,cv_H2O,mu,K,lamW,phi,g) = kante.Param
    #--

    #-- Zustandsvariablen
    mR = kante.y.mR
    eR = kante.y.eR
    #--

    (; KR,RR) = kante
    PR = KR.y.P
    TR = KR.y.T
    mnx = RR.y.m
    Pnx = RR.y.P
    Tnx = RR.y.T

    lam = lamda(mR,D,A,mu,K)

    dy[k] = -(mR^2-mnx^2)*2/(dx*Arho) - A*(PR-Pnx)*2/dx - lam/(2*D*Arho)*abs(mR)*mR - g*Arho*sin(phi); #-- mR
    dy[k+1] = eR -(0.5*cv_H2O*(abs(mR)*(Tnx-TR)+mR*(Tnx+TR)) + A/dx*2*lamW*(Tnx-TR)); #-- eR
end

function Kante!(dy,k,kante::mWRo_kante_mitte,t)
    #-- Parameter
    (; a2,A,D,dx,Arho,leit,K,mu,phi,g) = kante.Param
    #--

    #-- Zustandsvariablen
    P = kante.y.P
    m = kante.y.m
    T = kante.y.T
    #--

    lam = lamda(m,D,A,mu,K)

    (; RL,RR) = kante
    if typeof(kante.RL)==mWRo_kante_randL
        P_m12 = RL.KL.y.P; m_m12 = RL.y.mL; T_m12 = RL.KL.y.T
    else
        P_m12 = 0.5*(P + RL.y.P); m_m12 = 0.5*(m + RL.y.m); T_m12 = 0.5*(T + RL.y.T)
    end
    if typeof(kante.RR)==mWRo_kante_randR
        P_p12 = RR.KR.y.P;  m_p12 = RR.y.mR; T_p12 = RR.KR.y.T 
    else
        P_p12 = 0.5*(P + RR.y.P); m_p12 = 0.5*(m + RR.y.m); T_p12 = 0.5*(T + RR.y.T)
    end
    dy[k] = -a2/A*(m_p12-m_m12)/dx  #-- P
    dy[k+1] = -(m_p12^2-m_m12^2)/(dx*Arho) - A*(P_p12-P_m12)/dx - lam/(2*D*Arho)*abs(m)*m - g*Arho*sin(phi)  #-- m
    dy[k+2] = -1/Arho*m*ifxaorb(m,T-T_m12,T_p12-T)*2/dx + leit*2/(dx^2)*(T_m12-2*T+T_p12)   #-- T
end


function lamda(m,D,A,mu,K)
    Re = abs(m)*D/(A*mu); Re = max(Re,1.0e-6)
    lam = 0.25/(log10(K/(3.7*D)+5.74/exp(0.9*log(Re)))^2)
    return lam
end