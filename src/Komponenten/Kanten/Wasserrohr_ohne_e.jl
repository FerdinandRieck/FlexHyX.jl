Base.@kwdef mutable struct mWRoE_Param
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

Base.@kwdef mutable struct y_mWRoE
    Param::mWRoE_Param
    mL::Number = Param.m
    mR::Number = Param.m
    P
    _m
    T
end

Base.@kwdef mutable struct mWRoE_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWRoE_Param

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- Zustandsvariablen
    y = y_mWRoE(Param=Param,
               P = f(Vector(1:2:2*Param.nx), KL.y.P, KR.y.P, Param.nx), 
               T = f(Vector(1:2:2*Param.nx), KL.y.T, KR.y.T, Param.nx),
               _m = fill(Param.m,Param.nx)
               )

    #-- M-Matrix
    M::Array{Int} = [1; 1; ones(Int,3*Param.nx)]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWRoE_kante,t)
    #-- Parameter
    (; nx,dx,a2,leit,Arho,rho0,A,Di,cv_H2O,mu,K,lamW,phi,g,WENO,fluxPL,fluxPR,fluxmL,fluxmR,fluxTL,fluxTR) = kante.Param
    #--

    #-- Zustandsvariablen
    mL = kante.y.mL
    mR = kante.y.mR
    P = kante.y.P
    m = kante.y._m
    T = kante.y.T
    #--

    (; KL,KR,Z) = kante
    PL = KL.y.P
    TL = KL.y.T
    PR = KR.y.P
    TR = KR.y.T

    if WENO == true
        recover_weno!(P,fluxPL,fluxPR)
        recover_weno!(m,fluxmL,fluxmR)
        recover_weno!(T,fluxTL,fluxTR) #??? hier vieleicht nur recover(T), da bereits ubwind Diskertisierung mit ifxaorb(m[i],T[i]-fLT,fRT-T[i])??? 
    else
        recover!(P,fluxPL,fluxPR)
        recover!(m,fluxmL,fluxmR)
        recover!(T,fluxTL,fluxTR)
    end

    #-- Rohr links
    dy[k] = -(m[1]^2-mL^2)*2/(dx*Arho) - 1e5*A*(P[1]-PL)*2/dx - lambda(mL,Di,A,mu,K)/(2*Di*Arho)*abs(mL)*mL - g*Arho*sin(phi); #-- mL
    #-- Rohr rechts
    dy[k+1] = -(mR^2-m[end]^2)*2/(dx*Arho) - 1e5*A*(PR-P[end])*2/dx - lambda(mR,Di,A,mu,K)/(2*Di*Arho)*abs(mR)*mR - g*Arho*sin(phi); #-- mR
    
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
        dy[k+i+1] = -a2/A*(fRm-fLm)/dx*1e-5
        dy[k+i+1+nx] = -(fRm^2-fLm^2)/(dx*Arho) - 1e5*A*(fRP-fLP)/dx - lambda(m[i],Di,A,mu,K)/(2*Di*Arho)*abs(m[i])*m[i] - g*Arho*sin(phi)
        dy[k+i+1+nx*2] = -1/Arho*m[i]*ifxaorb(m[i],T[i]-fLT,fRT-T[i])*2/dx + leit*2/(dx^2)*(fLT-2*T[i]+fRT)
        if haskey(Z,"kA")
            if isa(Z["kA"],Number) 
                dy[k+i+1+nx*2] = dy[k+i+1+nx*2] - 4*Z["kA"]/(rho0*Di*cv_H2O)*(T[i]-Z["T_aussen"]) 
            end
            if isa(Z["kA"],Function)
                dy[k+i+1+nx*2] = dy[k+i+1+nx*2] - 4*Z["kA"](kante)/(rho0*Di*cv_H2O)*(T[i]-Z["T_aussen"])
            end
        end
    end
    
    #-- Energieflüsse berechnen und Zusatzinfos schreiben
    TRL = T[1] - (T[1]-T[2])/dx * - 0.5*dx 
    TRR = T[nx-1] - (T[nx-1]-T[nx])/dx * 1.5*dx
    eL = 1e-6*(cv_H2O*0.5*(abs(mL)*(TL-TRL)+mL*(TL+TRL)) + A/dx*2*lamW*(TL-TRL)) #-- eL
    eR = 1e-6*(cv_H2O*0.5*(abs(mR)*(TRR-TR)+mR*(TRR+TR)) + A/dx*2*lamW*(TRR-TR))
    Z["eL"] = eL; Z["eR"] = eR
end

function lambda(m,D,A,mu,K)
    Re = abs(m)*D/(A*mu); Re = max(Re,1.0e-6)
    lam = 0.25/(log10(K/(3.7*D)+5.74/exp(0.9*log(Re)))^2)
    return lam
end

function recover!(y,yL,yR)
    if length(y) > 1
        yL[1] = 2*y[1]-y[2];  yL[2:end] = y
        yR[end]= 2*y[end]-y[end-1]; yR[1:end-1] = y
    end
    return nothing
end

function recover_weno!(y,yL,yR)
    #-- Zellenmittelwerte auf Zellgrenzen interpolieren
    #-- L = upwind, R = downwind, WENO 3. Ordnung
    n = length(y);
    ep = 1.0e-6; p = 0.6;
    if n > 2
        yL[1] = 11/6*y[1]-7/6*y[2]+y[3]/3; yR[1] = yL[1]; #-- Randwerte
        yL[2] = y[1]/3+5/6*y[2]-y[3]/6;
        yL[n+1] = 11/6*y[n]-7/6*y[n-1]+y[n-2]/3; yR[n+1] = yL[n+1];
        yR[n] = y[n]/3+5/6*y[n-1]-y[n-2]/6; 
        for i=2:n-1
            uL = y[i]-y[i-1]; uC = y[i+1]-2*y[i]+y[i-1]; uR = y[i+1]-y[i]; uCC = y[i+1]-y[i-1];
            ISL = uL^2; ISC = 13/3*uC^2 +0.25*uCC^2; ISR = uR^2; aL = 0.25*(1/(ep+ISL))^p; aC = 0.5*(1/(ep+ISC))^p;
            aR = 0.25*(1/(ep+ISR))^p;
            suma = max(aL+aC+aR,eps(1.0)); wL = aL/suma; wC = aC/suma; wR = aR/suma;
            yR[i] = (0.5*wL+5/12*wC)*y[i-1] + (0.5*wL+2/3*wC+1.5*wR)*y[i] + (-wC/12-0.5*wR)*y[i+1];
            yL[i+1] = (-0.5*wL-wC/12)*y[i-1] + (1.5*wL+2/3*wC+0.5*wR)*y[i] + (5/12*wC+0.5*wR)*y[i+1];
        end 
    end
    return nothing
end



