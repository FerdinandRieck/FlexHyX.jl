Base.@kwdef mutable struct mWRo2_Param
    nx = 1
    L = 1.0
    dx = L/max(nx,1/2)
    Di = 0.1
    Dm = 0.11
    Da = 0.17
    A = pi*(Di/2)^2
    rho0 = 1000.0
    lamW = 0.6
    lamRohr1 = 0.13
    lamRohr2 = 0.040
    alpha_a = 10.0 #-- Wärmeübergang bei einer Luftbewegung von 0.1 m/s (arbeitswissenschaftliche Empfehlung für Büroraum)
    cv_H2O = 4182.0; #-- noch faglich
    mu = 1.0e-3
    Arho = A*rho0
    leit = lamW/(rho0*cv_H2O)
    g = 9.81
    a = 1414
    a2 = a^2
    K = 1e-5 #-- Rauheit
    phi = 0.0 #-- Neigungswinkel
    m = 0.0
    WENO = true
    m_dot = 1.0
    fluxPL = Array{Number}(undef, nx+1)
    fluxPR = Array{Number}(undef, nx+1)
    fluxmL = Array{Number}(undef, nx+1)
    fluxmR = Array{Number}(undef, nx+1)
    fluxTL = Array{Number}(undef, nx+1)
    fluxTR = Array{Number}(undef, nx+1)
end

#-- Rohrränder ---------------------------------
Base.@kwdef mutable struct y_mWRo2
    mL::Number = 0.0
    eL::Number = 0.0
    P
    _m
    T
    mR::Number = 0.0
    eR::Number = 0.0
end

Base.@kwdef mutable struct mWRo2_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWRo2_Param

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- Zustandsvariablen
    y = y_mWRo2(P = f(Vector(1:2:2*Param.nx), KL.y.P, KR.y.P, Param.nx), 
                T = f(Vector(1:2:2*Param.nx), KL.y.T, KR.y.T, Param.nx),
                _m = fill(Param.m,Param.nx)
                )

    #-- M-Matrix
    M::Array{Int} = [1; 0; ones(Int,3*Param.nx); 1; 0]

    #-- zusätzliche Infos
    Z::Dict
end


function Kante!(dy,k,kante::mWRo2_kante,t)
    #-- Parameter
    (; nx,dx,a2,a,leit,Arho,rho0,A,Di,cv_H2O,mu,K,lamW,phi,g,WENO,fluxPL,fluxPR,fluxmL,fluxmR,fluxTL,fluxTR) = kante.Param
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

    (; KL,KR,Z) = kante
    PL = KL.y.P
    TL = KL.y.T
    PR = KR.y.P
    TR = KR.y.T

    if WENO == true
        fluxPL, fluxPR = recover_weno(P)
        fluxmL, fluxmR = recover_weno(m)
        fluxTL, fluxTR = recover_weno(T)
    else
        fluxPL, fluxPR = recover(P)
        fluxmL, fluxmR = recover(m)
        fluxTL, fluxTR = recover(T)
    end
    #-- Rohr links
    dy[k] = -(m[1]^2-mL^2)*2/(dx*Arho) - 1e5*A*(P[1]-PL)*2/dx - lambda(mL,Di,A,mu,K)/(2*Di*Arho)*abs(mL)*mL - g*Arho*sin(phi); #-- mL
    TRL = 0.5*(3*T[1] - T[2]) #-- Extrapolation der Temp. auf Rohreinlauf
    dy[k+1] = eL - 1e-6*(cv_H2O*0.5*(abs(mL)*(TL-TRL)+mL*(TL+TRL)) + A/dx*2*lamW*(TL-T[1])) #-- eL  
    #-- Rohr mitte
    fRP = PL; fRm = mL; fRT = TL #-- linke Randbedingung
    fCRP = 0.5*(fluxPR[1]-PL); fCRm = 0.5*(fluxmR[1]-mL)  #-- Flusskorrekturen
    fluxTL[1] = TL; fluxTR[nx+1] = TR;
    for i = 1:nx
        fLP = fRP; fLm = fRm; fLT = fRT
        fCLP = fCRP; fCLm = fCRm  #-- Flusskorrekturen
        if i == nx
            fRP = PR; fRm = mR ; fRT = TR
            fCRP = 0.5*(PR-fluxPL[nx+1]); fCRm = 0.5*(mR-fluxmL[nx+1]);   #-- Flusskorrekturen
        else
            fRP = 0.5*(fluxPL[i+1]+fluxPR[i+1]); fRm = 0.5*(fluxmL[i+1]+fluxmR[i+1]); fRT = 0.5*(fluxTL[i+1]+fluxTR[i+1])
            fCRP = 0.5*(fluxPR[i+1]-fluxPL[i+1]); fCRm = 0.5*(fluxmR[i+1]-fluxmL[i+1])  #-- Flusskorrekturen
        end
        dy[k+i+1] = -a2/A*(fRm-fLm)/dx * 1e-5
        dy[k+i+1+nx] = -(fRm^2-fLm^2)/(dx*Arho) - 1e5*A*(fRP-fLP)/dx - lambda(m[i],Di,A,mu,K)/(2*Di*Arho)*abs(m[i])*m[i] - g*Arho*sin(phi)
        if WENO  
            dy[k+i+1] = dy[k+i+1] + a*(fCRP-fCLP)/dx #-- dämpft zu stark
            dy[k+i+1+nx] = dy[k+i+1+nx] + a*(fCRm-fCLm)/dx
        end
        dy[k+i+1+nx*2] =  -m[i]/Arho*ifxaorb(m[i],fluxTL[i+1]-fluxTL[i],fluxTR[i+1]-fluxTR[i])/dx + leit*2/(dx^2)*(fLT-2*T[i]+fRT)
        if haskey(Z,"kA")
            if isa(Z["kA"],Number) 
                dy[k+i+1+nx*2] = dy[k+i+1+nx*2] - 4*Z["kA"]/(rho0*Di*cv_H2O)*(T[i]-Z["T_aussen"]) 
            end
            if isa(Z["kA"],Function)
                dy[k+i+1+nx*2] = dy[k+i+1+nx*2] - 4*Z["kA"](kante)/(rho0*Di*cv_H2O)*(T[i]-Z["T_aussen"])
            end
        end
    end
    #-- Rohr rechts
    dy[k+3*nx+2] = -(mR^2-m[end]^2)*2/(dx*Arho) - 1e5*A*(PR-P[end])*2/dx - lambda(mR,Di,A,mu,K)/(2*Di*Arho)*abs(mR)*mR - g*Arho*sin(phi); #-- mR
    TRR = 0.5*(3*T[nx] - T[nx-1]) #-- Extrapolation der Temp. auf Rohrauslauf
    dy[k+3*nx+3] = eR - 1e-6*(cv_H2O*0.5*(abs(mR)*(TRR-TR)+mR*(TRR+TR)) + A/dx*2*lamW*(T[nx]-TR)) #-- eR
end