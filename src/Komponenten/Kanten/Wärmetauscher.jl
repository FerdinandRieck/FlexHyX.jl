Base.@kwdef mutable struct mWTaR_Param
    nx = 1
    L = 1.0
    dx = L/max(nx,1/2)
    Dü = 0.02
    Da = 0.1
    D = Da
    A = pi/4*abs(Da^2-Dü^2)
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
    m = 0.0
    WENO = true
    Richtung = "gleich"
    fluxPL = Array{Number}(undef, nx+1)
    fluxPR = Array{Number}(undef, nx+1)
    fluxmL = Array{Number}(undef, nx+1)
    fluxmR = Array{Number}(undef, nx+1)
    fluxTL = Array{Number}(undef, nx+1)
    fluxTR = Array{Number}(undef, nx+1)
end

Base.@kwdef mutable struct y_mWTaR
    Param::mWTaR_Param
    mL::Number = Param.m
    eL::Number = 0.0
    P
    _m
    T
    mR::Number = Param.m
    eR::Number = 0.0
end

Base.@kwdef mutable struct mWTaR_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWTaR_Param

    #-- Knoten links und rechts
    KL::Knoten
    KR::Knoten

    #-- Zustandsvariablen
    y = y_mWTaR(Param=Param,
                P = f(Vector(1:2:2*Param.nx), KL.y.P, KR.y.P, Param.nx), 
                T = f(Vector(1:2:2*Param.nx), KL.y.T, KR.y.T, Param.nx),
                _m = fill(Param.m,Param.nx)
                )

    #-- Rohr aussen
    R = 0

    #-- M-Matrix
    M::Array{Int} = [1; 0; ones(Int,3*Param.nx); 1; 0]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWTaR_kante,t)
    #-- Parameter
    (; nx,dx,a2,leit,Arho,rho0,A,Dü,D,cv_H2O,mu,K,lamW,phi,g,kA,WENO,Richtung,fluxPL,fluxPR,fluxmL,fluxmR,fluxTL,fluxTR) = kante.Param
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

    (; KL,KR,R) = kante
    PL = KL.y.P
    TL = KL.y.T
    PR = KR.y.P
    TR = KR.y.T
    T_aussen = R.y.T

    if WENO == true
        recover_weno!(P,fluxPL,fluxPR)
        recover_weno!(m,fluxmL,fluxmR)
        recover_weno!(T,fluxTL,fluxTR) #??? hier vieleicht nur recover!(T), da bereits ubwind Diskertisierung mit ifxaorb(m[i],T[i]-fLT,fRT-T[i])??? 
    else
        recover!(P,fluxPL,fluxPR)
        recover!(m,fluxmL,fluxmR)
        recover!(T,fluxTL,fluxTR)
    end

    if (Richtung == "gegen") && nx > 1
        T_aussen = reverse(T_aussen)
    end

    #-- Rohr links
    dy[k] = -(m[1]^2-mL^2)*2/(dx*Arho) - A*(P[1]-PL)*2/dx - lambda(mL,D,A,mu,K)/(2*D*Arho)*abs(mL)*mL - g*Arho*sin(phi)  #-- mL
    TRL = T[1] - (T[1]-T[2])/dx * -0.5*dx
    dy[k+1] = eL -(0.5*cv_H2O*(abs(mL)*(TL-TRL)+mL*(TL+TRL)) + A/dx*2*lamW*(TL-T[1])); #-- eL
    
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
        dy[k+i+1+nx*2] = -1/Arho*m[i]*ifxaorb(m[i],T[i]-fLT,fRT-T[i])*2/dx + leit*2/(dx^2)*(fLT-2*T[i]+fRT) -  kA*pi*Dü/(Arho*cv_H2O)*(T[i]-T_aussen[i])  #-- T
    end

    #-- Rohr rechts
    dy[k+3*nx+2] = -(mR^2-m[end]^2)*2/(dx*Arho) - A*(PR-P[end])*2/dx - lambda(mR,D,A,mu,K)/(2*D*Arho)*abs(mR)*mR - g*Arho*sin(phi)  #-- mR
    TRR = T[nx-1] - (T[nx-1]-T[nx])/dx * 1.5*dx
    dy[k+3*nx+3] = eR -(0.5*cv_H2O*(abs(mR)*(TRR-TR)+mR*(TRR+TR)) + A/dx*2*lamW*(T[end]-TR)) #-- eR
end