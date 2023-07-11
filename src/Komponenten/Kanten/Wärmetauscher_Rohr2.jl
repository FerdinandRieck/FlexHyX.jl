Base.@kwdef mutable struct mWTaR2_Param
    nx = 1
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
    PL = 10
    PR = 0
    TL = 10
    TR = 0
    P_vec = PL:(PR-PL)*0.5/nx:PR
    T_vec = TL:(TR-TL)*0.5/nx:TR
    P = Vector(P_vec[2:2:end-1])
    T = Vector(T_vec[2:2:end-1])
    kA = 380.0
end

#-- Rohrränder ---------------------------------
Base.@kwdef mutable struct y_mWTaR2
    Param::mWTaR2_Param
    mL::Number = 0.0
    eL::Number = 0.0
    P = Param.P
    _m = zeros(Param.nx)
    T = Param.T
    mR::Number = 0.0
    eR::Number = 0.0
end

Base.@kwdef mutable struct mWTaR2_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWTaR2_Param

    #-- Zustandsvariablen
    y = y_mWTaR2(Param=Param)

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- Rohr aussen
    RA = 0

    #-- M-Matrix
    M::Array{Int} = [1; 0; ones(Int,3*Param.nx); 1; 0]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWTaR2_kante,t)
    #-- Parameter
    (; nx,dx,a2,leit,Arho,A,D,cv_H2O,mu,K,lamW,phi,g,kA) = kante.Param
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

    (; KL,KR,RA,Z) = kante
    PL = KL.y.P
    TL = KL.y.T
    PR = KR.y.P
    TR = KR.y.T
    T_aussen = RA.y.T

    if haskey(Z,"WENO") 
        if Z["WENO"] == true
            fluxPL, fluxPR = recover_weno(P)
            fluxmL, fluxmR = recover_weno(m)
            fluxTL, fluxTR = recover_weno(T)
        end
        if Z["WENO"] == false
            fluxPL, fluxPR = recover(P)
            fluxmL, fluxmR = recover(m)
            fluxTL, fluxTR = recover(T)
        end
    end

    #-- Rohr links
    dy[k] = -(m[1]^2-mL^2)*2/(dx*Arho) - A*(P[1]-PL)*2/dx - lamda(mL,D,A,mu,K)/(2*D*Arho)*abs(mL)*mL - g*Arho*sin(phi)  #-- mL
    dy[k+1] = eL -(0.5*cv_H2O*(abs(mL)*(TL-T[1])+mL*(TL+T[1])) + A/dx*2*lamW*(TL-T[1]))  #-- eL
    
    #-- Rohr mitte
    fRP = PL; fRm = mL; fRT = TL #-- linke Randbedingung
    for i = 1:nx
        fLP = fRP; fLm = fRm; fLT = fRT
        if i == nx
            fRP = PR; fRm = mR ; fRT = TR
        else
            fRP = 0.5*(fluxPL[i+1]+fluxPR[i+1]) 
            fRm = 0.5*(fluxmL[i+1]+fluxmR[i+1]) 
            fRT = 0.5*(fluxTL[i+1]+fluxTR[i+1])
        end
        dy[k+i+1] = -a2/A*(fRm-fLm)/dx  #-- P 
        dy[k+i+1+nx] = -(fRm^2-fLm^2)/(dx*Arho) - A*(fRP-fLP)/dx - lamda(m[i],D,A,mu,K)/(2*D*Arho)*abs(m[i])*m[i] - g*Arho*sin(phi)  #-- m #!!! Nicht klar ob Winkel Angegeben werden darf, wegen GL.88 (Herleitung aus Text Rohr (Wasser))
        dy[k+i+1+nx*2] = -1/Arho*m[i]*ifxaorb(m[i],T[i]-fLT,fRT-T[i])*2/dx + leit*2/(dx^2)*(fLT-2*T[i]+fRT) - kA/(Arho*D)*(T[i]-T_aussen[i])  #-- T
    end

    #-- Rohr rechts
    dy[k+3*nx+2] = -(mR^2-m[end]^2)*2/(dx*Arho) - A*(PR-P[end])*2/dx - lamda(mR,D,A,mu,K)/(2*D*Arho)*abs(mR)*mR - g*Arho*sin(phi)  #-- mR
    dy[k+3*nx+3] = eR -(0.5*cv_H2O*(abs(mR)*(T[end]-TR)+mR*(T[end]+TR)) + A/dx*2*lamW*(T[end]-TR)) #-- eR
end

function lamda(m,D,A,mu,K)
    Re = abs(m)*D/(A*mu); Re = max(Re,1.0e-6)
    lam = 0.25/(log10(K/(3.7*D)+5.74/exp(0.9*log(Re)))^2)
    return lam
end

function recover(y)
    yL = [2*y[1]-y[2]; y]; yR = [y; 2*y[end]-y[end-1]] 
    return yL, yR
end

function recover_weno(y)
    #-- Zellenmittelwerte auf Zellgrenzen interpolieren
    #-- L = upwind, R = downwind, WENO 3. Ordnung
    n = length(y);
    yL = Array{Number}(undef, n+1); yR = Array{Number}(undef, n+1) 
    yL[1] = 11/6*y[1]-7/6*y[2]+y[3]/3; yR[1] = yL[1]; #-- Randwerte
    yL[2] = y[1]/3+5/6*y[2]-y[3]/6;
    yL[n+1] = 11/6*y[n]-7/6*y[n-1]+y[n-2]/3; yR[n+1] = yL[n+1];
    yR[n] = y[n]/3+5/6*y[n-1]-y[n-2]/6; 
    for i=2:n-1
        yR[i], yL[i+1] = weno3(y[i-1:i+1]); 
    end
    return yL, yR 
end

function weno3(y) #-- y = [y1,y2,y3]
    ep = 1.0e-6; p = 0.6;
    uL = y[2]-y[1]; uC = y[3]-2*y[2]+y[1]; uR = y[3]-y[2]; uCC = y[3]-y[1];
    ISL = uL^2; ISC = 13/3*uC^2 +0.25*uCC^2; ISR = uR^2; aL = 0.25*(1/(ep+ISL))^p; aC = 0.5*(1/(ep+ISC))^p;
    aR = 0.25*(1/(ep+ISR))^p;
    suma = max(aL+aC+aR,eps(1.0)); wL = aL/suma; wC = aC/suma; wR = aR/suma;
    y12 = (0.5*wL+5/12*wC)*y[1] + (0.5*wL+2/3*wC+1.5*wR)*y[2] + (-wC/12-0.5*wR)*y[3];
    y23 = (-0.5*wL-wC/12)*y[1] + (1.5*wL+2/3*wC+0.5*wR)*y[2] + (5/12*wC+0.5*wR)*y[3];
    return y12, y23
end