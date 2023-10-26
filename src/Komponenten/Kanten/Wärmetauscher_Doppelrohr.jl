Base.@kwdef mutable struct mWTaR_Param
    nx = 1
    L = 1.0
    dx = L/max(nx,1/2)
    rho0 = 1000
    Di = 0.02
    Dm = 0.022
    Da = 0.1
    Aa = pi/4*(Da^2-Dm^2)
    Ai = pi/4*Di^2
    Aarho = Aa*rho0
    Airho = Ai*rho0
    lamW = 0.6
    lamRohr = 401.0 #-- Wärmleitung Kupfer
    kA = 2000.0
    cv_H2O = 4182.0; #-- noch faglich
    leit = lamW/(rho0*cv_H2O)
    mu = 1.0e-3
    WENO = true
    Richtung = "gleich"
    Ringspalt = true
    g = 9.81
    a = 1414
    a2 = a^2
    K = 1e-5 #-- Rauheit
    phi = 0.0 #-- Neigungswinkel
    m = 0.0
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

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

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
    (; nx,dx,a2,leit,Ai,Aa,Aarho,Airho,Di,Dm,Da,cv_H2O,mu,K,lamW,phi,g,WENO,Richtung,fluxPL,fluxPR,fluxmL,fluxmR,fluxTL,fluxTR,Ringspalt) = kante.Param
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

    (; KL,KR,R,Z) = kante
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

    if Ringspalt==true
        if isa(Z["kA"],Function)
            kA = Z["kA"](kante)
        else
            kA = Z["kA"]
        end
        Arho = Aarho
        A = Aa
        D = Da
    elseif isa(Ringspalt,Int)
        kA = Z["kA"]
        Arho = Airho
        A = Ai
        D = Di
    end

    #-- Rohr links
    dy[k] = -(m[1]^2-mL^2)*2/(dx*Arho) - 1e5*A*(P[1]-PL)*2/dx - lambda(mL,D,A,mu,K)/(2*D*Arho)*abs(mL)*mL - g*Arho*sin(phi)  #-- mL
    TRL = T[1] - (T[1]-T[2])/dx * -0.5*dx
    dy[k+1] = eL - 1e-6*(cv_H2O*0.5*(abs(mL)*(TL-TRL)+mL*(TL+TRL)) + A/dx*2*lamW*(TL-TRL)) #-- eL 
    
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
        dy[k+i+1] = -a2/A*(fRm-fLm)/dx*1e-5  #-- P 
        dy[k+i+1+nx] = -(fRm^2-fLm^2)/(dx*Arho) - 1e5*A*(fRP-fLP)/dx - lambda(m[i],D,A,mu,K)/(2*D*Arho)*abs(m[i])*m[i] - g*Arho*sin(phi)  #-- m 
        dy[k+i+1+nx*2] = -1/Arho*m[i]*ifxaorb(m[i],T[i]-fLT,fRT-T[i])*2/dx + leit*2/(dx^2)*(fLT-2*T[i]+fRT) -  kA*pi*Dm/(Arho*cv_H2O)*(T[i]-T_aussen[i])  #-- T
    end

    #-- Rohr rechts
    dy[k+3*nx+2] = -(mR^2-m[end]^2)*2/(dx*Arho) - 1e5*A*(PR-P[end])*2/dx - lambda(mR,D,A,mu,K)/(2*D*Arho)*abs(mR)*mR - g*Arho*sin(phi)  #-- mR
    TRR = T[nx-1] - (T[nx-1]-T[nx])/dx * 1.5*dx
    dy[k+3*nx+3] = eR - 1e-6*(cv_H2O*0.5*(abs(mR)*(TRR-TR)+mR*(TRR+TR)) + A/dx*2*lamW*(TRR-TR)) #-- eR
end