Base.@kwdef mutable struct mWTaR_Param
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
    P_vec = PL:(PR-PL)*0.5/nx:PR #[fill(PL,nx);fill(PR,nx+1)]
    T_vec = [fill(TL,nx);fill(TR,nx+1)] #TL:(TR-TL)*0.5/nx:TR
    P = Vector(P_vec[2:2:end-1])
    T = Vector(T_vec[2:2:end-1])
    kA = 380.0
    WENO = true
    Richtung = "gegen"
end

#-- Rohrränder ---------------------------------
Base.@kwdef mutable struct y_mWTaR
    Param::mWTaR_Param
    mL::Number = 0.0
    eL::Number = 0.0
    P = Param.P
    _m = zeros(Param.nx)
    T = Param.T
    mR::Number = 0.0
    eR::Number = 0.0
end

Base.@kwdef mutable struct mWTaR_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWTaR_Param

    #-- Zustandsvariablen
    y = y_mWTaR(Param=Param)

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

function Kante!(dy,k,kante::mWTaR_kante,t)
    #-- Parameter
    (; nx,dx,a2,leit,Arho,A,D,cv_H2O,mu,K,lamW,phi,g,kA,WENO,Richtung) = kante.Param
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
        fluxPL, fluxPR = recover_weno(P)
        fluxmL, fluxmR = recover_weno(m)
        fluxTL, fluxTR = recover_weno(T)
    else
        fluxPL, fluxPR = recover(P)
        fluxmL, fluxmR = recover(m)
        fluxTL, fluxTR = recover(T)
    end

    if Richtung == "gegen"
        T_aussen = reverse(T_aussen)
    end


    #-- Rohr links
    dy[k] = -(m[1]^2-mL^2)*2/(dx*Arho) - A*(P[1]-PL)*2/dx - lambda(mL,D,A,mu,K)/(2*D*Arho)*abs(mL)*mL - g*Arho*sin(phi)  #-- mL
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
        dy[k+i+1+nx] = -(fRm^2-fLm^2)/(dx*Arho) - A*(fRP-fLP)/dx - lambda(m[i],D,A,mu,K)/(2*D*Arho)*abs(m[i])*m[i] - g*Arho*sin(phi)  #-- m #!!! Nicht klar ob Winkel Angegeben werden darf, wegen GL.88 (Herleitung aus Text Rohr (Wasser))
        dy[k+i+1+nx*2] = -1/Arho*m[i]*ifxaorb(m[i],T[i]-fLT,fRT-T[i])*2/dx + leit*2/(dx^2)*(fLT-2*T[i]+fRT) - kA/(Arho*D)*(T[i]-T_aussen[i])  #-- T
    end

    #-- Rohr rechts
    dy[k+3*nx+2] = -(mR^2-m[end]^2)*2/(dx*Arho) - A*(PR-P[end])*2/dx - lambda(mR,D,A,mu,K)/(2*D*Arho)*abs(mR)*mR - g*Arho*sin(phi)  #-- mR
    dy[k+3*nx+3] = eR -(0.5*cv_H2O*(abs(mR)*(T[end]-TR)+mR*(T[end]+TR)) + A/dx*2*lamW*(T[end]-TR)) #-- eR
end