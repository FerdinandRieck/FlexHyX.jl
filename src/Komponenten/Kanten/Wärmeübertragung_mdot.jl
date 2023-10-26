Base.@kwdef mutable struct mWTM_Param
    nx = 1
    L = 1.0
    dx = L/max(nx,1/2)
    rho0 = 1000.0
    Di = 0.02
    Dm = 0.022
    Da = 0.026
    A = pi/4*Di^2
    lamW = 0.6
    lamRohr1 = 160.0 #-- Wärmleitung Aluminium
    lamRohr2 = 0.13 #-- Wärmleitung Holz
    alpha_a = 10.0 #-- Wärmeübergang bei einer Luftbewegung von 0.1 m/s (arbeitswissenschaftliche Empfehlung für Büroraum)
    cv_H2O = 4182.0 #-- noch faglich
    Arho = A*rho0
    leit = lamW/(rho0*cv_H2O)
    mu = 1.0e-3
    kA = 100.0
    WENO = true
    fluxTL = Array{Number}(undef, nx+1)
    fluxTR = Array{Number}(undef, nx+1)
    m_dot = 1.0
end

Base.@kwdef mutable struct y_mWTM
    m::Number = 0.0
    eL::Number = 0.0
    eR::Number = 0.0
    T
end


Base.@kwdef mutable struct mWTM_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWTM_Param

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten
    
    #-- Knoten aussen
    KA = []

    #-- Zustandsvariablen
    y = y_mWTM(T = f(Vector(1:2:2*Param.nx), KL.y.T, KR.y.T, Param.nx))

    #-- M-Matrix
    M::Array{Int} = [0; 0; 0; ones(Int,Param.nx)]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWTM_kante,t)
    #-- Parameter
    (; nx,dx,leit,A,Arho,cv_H2O,Di,lamW,WENO,fluxTL,fluxTR,rho0) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    eL = kante.y.eL
    eR = kante.y.eR
    T = kante.y.T
    #--

    (; KL,KR,KA,Z) = kante
    TL = KL.y.T
    TR = KR.y.T
    T_aussen = KA.y.T
    mL = KL.in[1].y.m

    if WENO == true
        recover_weno!(T,fluxTL,fluxTR) 
    else
        recover!(T,fluxTL,fluxTR)
    end

    io = 1.0; if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end    

    #-- Rohr links
    dy[k] = m - mL #-- m    
    TRL = T[1] - (T[1]-T[2])/dx * -0.5*dx
    dy[k+1] = eL - 1e-6*(0.5*cv_H2O*(abs(m)*(TL-TRL)+m*(TL+TRL)) + A/dx*2*lamW*(TL-T[1])); #-- eL
    #-- Rohr rechts
    TRR = T[nx-1] - (T[nx-1]-T[nx])/dx * 1.5*dx
    dy[k+2] = eR - 1e-6*(0.5*cv_H2O*(abs(m)*(TRR-TR)+m*(TRR+TR)) + A/dx*2*lamW*(T[end]-TR)) #-- eR

    #-- Rohr mitte 
    fRT = TL #-- linke Randbedingung
    for i = 1:nx
        fLT = fRT
        if i == nx
            fRT = TR
        else 
            fRT = 0.5*(fluxTL[i+1]+fluxTR[i+1])
        end
        dy[k+i+2] = -1/Arho*m*ifxaorb(m,T[i]-fLT,fRT-T[i])*2/dx + leit*2/(dx^2)*(fLT-2*T[i]+fRT)  #-- T
        if haskey(Z,"kA")
            if isa(Z["kA"],Number) 
                dy[k+i+2] = dy[k+i+2] - 4*Z["kA"]/(rho0*Di*cv_H2O)*(T[i]-T_aussen) 
            end
            if isa(Z["kA"],Function)
                dy[k+i+2] = dy[k+i+2] - 4*Z["kA"](kante)/(rho0*Di*cv_H2O)*(T[i]-T_aussen)
            end
        end
    end
end