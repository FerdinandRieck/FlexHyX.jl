Base.@kwdef mutable struct mWRoD_Param
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
    m = 0.0
    WENO = true
    fluxPL = Array{Number}(undef, nx+1)
    fluxPR = Array{Number}(undef, nx+1)
    fluxmL = Array{Number}(undef, nx+1)
    fluxmR = Array{Number}(undef, nx+1)
    fluxTL = Array{Number}(undef, nx+1)
    fluxTR = Array{Number}(undef, nx+1)
end

Base.@kwdef mutable struct y_mWRoD
    Param::mWRoD_Param
    mL::Number = Param.m
    P
    _m
    mR::Number = Param.m
end

Base.@kwdef mutable struct mWRoD_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWRoD_Param

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- Zustandsvariablen
    y = y_mWRoD(Param=Param,
               P = f(Vector(1:2:2*Param.nx), KL.y.P, KR.y.P, Param.nx), 
               _m = fill(Param.m,Param.nx)
               )

    #-- M-Matrix
    M::Array{Int} = [1; ones(Int,2*Param.nx); 1]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWRoD_kante,t)
    #-- Parameter
    (; nx,dx,a2,leit,Arho,rho0,A,D,cv_H2O,mu,K,lamW,phi,g,WENO,fluxPL,fluxPR,fluxmL,fluxmR,fluxTL,fluxTR) = kante.Param
    #--

    #-- Zustandsvariablen
    mL = kante.y.mL
    P = kante.y.P
    m = kante.y._m
    mR = kante.y.mR
    #--

    (; KL,KR) = kante
    PL = KL.y.P
    PR = KR.y.P

    if WENO == true
        recover_weno!(P,fluxPL,fluxPR)
        recover_weno!(m,fluxmL,fluxmR)
    else
        recover!(P,fluxPL,fluxPR)
        recover!(m,fluxmL,fluxmR)
    end

    #-- Rohr links
    dy[k] = -(m[1]^2-mL^2)*2/(dx*Arho) - A*(P[1]-PL)*2/dx - lambda(mL,D,A,mu,K)/(2*D*Arho)*abs(mL)*mL - g*Arho*sin(phi); #-- mL

    #-- Rohr mitte
    fRP = PL; fRm = mL; #-- linke Randbedingung
    for i = 1:nx
        fLP = fRP; fLm = fRm; 
        if i == nx
            fRP = PR; fRm = mR ; 
        else
            fRP = 0.5*(fluxPL[i+1]+fluxPR[i+1]) 
            fRm = 0.5*(fluxmL[i+1]+fluxmR[i+1]) 
        end
        dy[k+i] = -a2/A*(fRm-fLm)/dx
        dy[k+i+nx] = -(fRm^2-fLm^2)/(dx*Arho) - A*(fRP-fLP)/dx - lambda(m[i],D,A,mu,K)/(2*D*Arho)*abs(m[i])*m[i] - g*Arho*sin(phi)
    end
    
    #-- Rohr rechts
    dy[k+2*nx+1] = -(mR^2-m[end]^2)*2/(dx*Arho) - A*(PR-P[end])*2/dx - lambda(mR,D,A,mu,K)/(2*D*Arho)*abs(mR)*mR - g*Arho*sin(phi); #-- mR
end
