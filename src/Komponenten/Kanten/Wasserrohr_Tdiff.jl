Base.@kwdef mutable struct mWRoTd_Param
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
Base.@kwdef mutable struct y_mWRoTd
    mL::Number = 0.0
    eL::Number = 0.0
    P
    _m
    T
    mR::Number = 0.0
    eR::Number = 0.0
end

function f(x,L,R,nx)
    y = L .+ x*0.5*(R-L)/nx
end

Base.@kwdef mutable struct mWRoTd_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWRoTd_Param

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- Zustandsvariablen
    y = y_mWRoTd(P = f(Vector(1:2:2*Param.nx), KL.y.P, KR.y.P, Param.nx), 
                T = f(Vector(1:2:2*Param.nx), KL.y.T, KR.y.T, Param.nx),
                _m = fill(Param.m,Param.nx)
                )

    #-- M-Matrix
    M::Array{Int} = [1; 0; ones(Int,3*Param.nx); 1; 0]

    #-- zusätzliche Infos
    Z::Dict
end


function Kante!(dy,k,kante::mWRoTd_kante,t)
    #-- Parameter
    (; nx,dx,a2,leit,Arho,rho0,A,Di,cv_H2O,mu,K,lamW,phi,g,WENO,fluxPL,fluxPR,fluxmL,fluxmR,fluxTL,fluxTR) = kante.Param
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
    dy[k] = -(m[1]^2-mL^2)*2/(dx*Arho) - 1e5*A*(P[1]-PL)*2/dx - lamda(mL,Di,A,mu,K)/(2*Di*Arho)*abs(mL)*mL - g*Arho*sin(phi); #-- mL
    TRL = T[1] - (T[1]-T[2])/dx * - 0.5*dx
    dy[k+1] =  eL - cv_H2O*mL*(TL-TRL) - A/dx*2*lamW*(TL-TRL) 
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
        dy[k+i+1+nx] = -(fRm^2-fLm^2)/(dx*Arho) - 1e5*A*(fRP-fLP)/dx - lamda(m[i],Di,A,mu,K)/(2*Di*Arho)*abs(m[i])*m[i] - g*Arho*sin(phi)
        if WENO  
            #dy[k+i+1] = dy[k+i+1] + a*(fCRP-fCLP)/dx #-- dämpft zu stark
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
    dy[k+3*nx+2] = -(mR^2-m[end]^2)*2/(dx*Arho) - 1e5*A*(PR-P[end])*2/dx - lamda(mR,Di,A,mu,K)/(2*Di*Arho)*abs(mR)*mR - g*Arho*sin(phi); #-- mR
    TRR = T[nx-1] - (T[nx-1]-T[nx])/dx * 1.5*dx
    dy[k+3*nx+3] = eR - cv_H2O*mR*(TRR-TR) - A/dx*2*lamW*(TRR-TR) 
end

function lamda(m,D,A,mu,K)
    Re = abs(m)*D/(A*mu); Re = max(Re,1.0e-6)
    lam = 0.25/(log10(K/(3.7*D)+5.74/exp(0.9*log(Re)))^2)
    #println(lam)
    #lam = 0.00001
    return lam
end

function recover(y)
    if length(y) > 1
        yL = [2*y[1]-y[2]; y]; yR = [y; 2*y[end]-y[end-1]] 
    else
        yL = [y[1]; y[1]]; yR = [y[end]; y[end]]
    end
    return yL, yR
end

function recover_weno(y)
    #-- Zellenmittelwerte auf Zellgrenzen interpolieren
    #-- L = upwind, R = downwind, WENO 3. Ordnung
    n = length(y);
    yL = Array{Number}(undef, n+1); yR = Array{Number}(undef, n+1) 
    #yL[1] = 11/6*y[1]-7/6*y[2]+y[3]/3; yR[1] = yL[1]; #-- Randwerte
    #yL[2] = y[1]/3+5/6*y[2]-y[3]/6;
    #yL[n+1] = 11/6*y[n]-7/6*y[n-1]+y[n-2]/3; yR[n+1] = yL[n+1];
    #yR[n] = y[n]/3+5/6*y[n-1]-y[n-2]/6; 
    yL[1] = 2*y[1]-y[2]; yR[1] = yL[1]; yL[2] = y[1]
    yR[n+1] = 2*y[n]-y[n-1]; yL[n+1] = yR[n+1]; yR[n] = y[n] 
    for i=2:n-1
        yR[i], yL[i+1] = weno3(y[i-1:i+1]); 
    end
    return yL, yR 
end

function weno3(y) #-- y = [y1,y2,y3]
    ep = 1.0e-6; p = 0.6;
    uL = y[2]-y[1]; uC = y[3]-2*y[2]+y[1]; uR = y[3]-y[2]; uCC = y[3]-y[1];
    ISL = uL^2; ISC = 13/3*uC^2 +0.25*uCC^2; ISR = uR^2; 
    aL = 0.25*(1/(ep+ISL))^p; aC = 0.5*(1/(ep+ISC))^p;  aR = 0.25*(1/(ep+ISR))^p;
    suma = max(aL+aC+aR,eps(1.0)); wL = aL/suma; wC = aC/suma; wR = aR/suma;
    y12 = (0.5*wL+5/12*wC)*y[1] + (0.5*wL+2/3*wC+1.5*wR)*y[2] + (-wC/12-0.5*wR)*y[3];
    y23 = (-0.5*wL-wC/12)*y[1] + (1.5*wL+2/3*wC+0.5*wR)*y[2] + (5/12*wC+0.5*wR)*y[3];
    return y12, y23
end
