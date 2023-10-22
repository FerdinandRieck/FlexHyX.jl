Base.@kwdef mutable struct mWRoM_Param
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
    WENO = true
    fluxTL = Array{Number}(undef, nx+1)
    fluxTR = Array{Number}(undef, nx+1)
    m_dot = 1.0
end

Base.@kwdef mutable struct y_mWRoM
    m::Number = 0.0
    eL::Number = 0.0
    eR::Number = 0.0
    T
end

Base.@kwdef mutable struct mWRoM_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWRoM_Param

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- Zustandsvariablen
    y = y_mWRoM(T = f(Vector(1:2:2*Param.nx), KL.y.T, KR.y.T, Param.nx))

    #-- M-Matrix
    M::Array{Int} = [0; 0; 0; ones(Int,Param.nx)]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWRoM_kante,t)
    #-- Parameter
    (; nx,dx,leit,A,Arho,rho0,Di,cv_H2O,lamW,WENO,fluxTL,fluxTR) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    eL = kante.y.eL
    eR = kante.y.eR
    T = kante.y.T
    #--

    (; KL,KR,Z) = kante
    TL = KL.y.T
    TR = KR.y.T
    mL = KL.in[1].y.m

    if WENO == true
        recover_weno!(T,fluxTL,fluxTR) #??? hier vieleicht nur recover(T), da bereits ubwind Diskertisierung mit ifxaorb(m[i],T[i]-fLT,fRT-T[i])??? 
    else
        recover!(T,fluxTL,fluxTR)
    end

    #-- Rohr links
    dy[k] = m - mL #-- m
    TRL = T[1] - (T[1]-T[2])/dx * - 0.5*dx
    dy[k+1] = eL - 1e-6*(0.5*cv_H2O*(abs(m)*(TL-TRL)+m*(TL+TRL)) + A/dx*2*lamW*(TL-TRL)) #-- eL
    #-- Rohr rechts
    TRR = T[nx-1] - (T[nx-1]-T[nx])/dx * 1.5*dx
    dy[k+2] = eR - 1e-6*(0.5*cv_H2O*(abs(m)*(TRR-TR)+m*(TRR+TR)) + A/dx*2*lamW*(TRR-TR)) #-- eR

    #-- Rohr mitte
    fRT = TL #-- linke Randbedingung
    for i = 1:nx
        fLT = fRT
        if i == nx
            fRT = TR
        else
            fRT = 0.5*(fluxTL[i+1]+fluxTR[i+1])
        end
        dy[k+i+2] = -1/Arho*m*ifxaorb(m,T[i]-fLT,fRT-T[i])*2/dx + leit*2/(dx^2)*(fLT-2*T[i]+fRT)
        if haskey(Z,"kA")
            if isa(Z["kA"],Number) 
                dy[k+i+2] = dy[k+i+2] - 4*Z["kA"]/(rho0*Di*cv_H2O)*(T[i]-Z["T_aussen"]) 
            end
            if isa(Z["kA"],Function)
                dy[k+i+2] = dy[k+i+2] - 4*Z["kA"](kante)/(rho0*Di*cv_H2O)*(T[i]-Z["T_aussen"])
            end
        end
    end
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
    if n > 1
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