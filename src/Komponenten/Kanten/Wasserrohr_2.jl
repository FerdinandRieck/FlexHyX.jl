Base.@kwdef mutable struct mWRo2_Param
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
    P_vec = PL:(PR-PL)*0.5/nx:PR # [fill(PL,nx);fill(PR,nx+1)]
    T_vec = TL:(TR-TL)*0.5/nx:TR
    P = Vector(P_vec[2:2:end-1])
    T = Vector(T_vec[2:2:end-1])
end

#-- Rohrränder ---------------------------------
Base.@kwdef mutable struct y_mWRo2
    Param::mWRo2_Param
    mL::Number = 0.0
    eL::Number = 0.0
    P = Param.P
    _m = zeros(Param.nx)
    T = Param.T
    mR::Number = 0.0
    eR::Number = 0.0
end

Base.@kwdef mutable struct mWRo2_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWRo2_Param

    #-- Zustandsvariablen
    y = y_mWRo2(Param=Param)

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- M-Matrix
    M::Array{Int} = [1; 0; ones(Int,3*Param.nx); 1; 0]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWRo2_kante,t)
    #-- Parameter
    #(; nx,dx,a2,leit,Arho,A,D,cv_H2O,mu,K,lamW,phi,g) = kante.Param
    #--
#=
    #-- Zustandsvariablen
    mL = kante.y.mL
    eL = kante.y.eL
    P = kante.y.P
    m = kante.y._m
    T = kante.y.T
    mR = kante.y.mR
    eR = kante.y.eR
    #--

    (; KL,KR) = kante
    PL = KL.y.P
    TL = KL.y.T
    PR = KR.y.P
    TR = KR.y.T
=#

    #-- Rohr links
    dy[k] = -(kante.y._m[1]^2-kante.y.mL^2)*2/(kante.Param.dx*kante.Param.Arho) - kante.Param.A*(kante.y.P[1]-kante.KL.y.P)*2/kante.Param.dx - lamda(kante.y.mL,kante.Param.D,kante.Param.A,kante.Param.mu,kante.Param.K)/(2*kante.Param.D*kante.Param.Arho)*abs(kante.y.mL)*kante.y.mL - kante.Param.g*kante.Param.Arho*sin(kante.Param.phi); #-- mL
    dy[k+1] = kante.y.eL -(0.5*kante.Param.cv_H2O*(abs(kante.y.mL)*(kante.KL.y.T-kante.y.T[1])+kante.y.mL*(kante.KL.y.T+kante.y.T[1])) + kante.Param.A/kante.Param.dx*2*kante.Param.lamW*(kante.KL.y.T-kante.y.T[1])); #-- eL
    
    #-- Rohr mitte
    for i = 1:kante.Param.nx
        if i==1
            P_m12 = kante.KL.y.P; m_m12 = kante.y.mL; T_m12 = kante.KL.y.T;
        else
            P_m12 = 0.5*(kante.y.P[i]+kante.y.P[i-1]); m_m12 = 0.5*(kante.y._m[i]+kante.y._m[i-1]); T_m12 = 0.5*(kante.y.T[i]+kante.y.T[i-1]);
        end
        if i==kante.Param.nx
            P_p12 = kante.KR.y.P; m_p12 = kante.y.mR; T_p12 = kante.KR.y.T;
        else
            P_p12 = 0.5*(kante.y.P[i]+kante.y.P[i+1]); m_p12 = 0.5*(kante.y._m[i]+kante.y._m[i+1]); T_p12 = 0.5*(kante.y.T[i]+kante.y.T[i+1]);
        end
        dy[k+i+1] = -kante.Param.a2/kante.Param.A*(m_p12-m_m12)/kante.Param.dx
        dy[k+i+1+kante.Param.nx] = -(m_p12^2-m_m12^2)/(kante.Param.dx*kante.Param.Arho) - kante.Param.A*(P_p12-P_m12)/kante.Param.dx - lamda(kante.y._m[i],kante.Param.D,kante.Param.A,kante.Param.mu,kante.Param.K)/(2*kante.Param.D*kante.Param.Arho)*abs(kante.y._m[i])*kante.y._m[i] - kante.Param.g*kante.Param.Arho*sin(kante.Param.phi)
        dy[k+i+1+kante.Param.nx*2] = -1/kante.Param.Arho*kante.y._m[i]*ifxaorb(kante.y._m[i],kante.y.T[i]-T_m12,T_p12-kante.y.T[i])*2/kante.Param.dx + kante.Param.leit*2/(kante.Param.dx^2)*(T_m12-2*kante.y.T[i]+T_p12)
    end

    #-- Rohr rechts
    dy[k+3*kante.Param.nx+2] = -(kante.y.mR^2-kante.y._m[end]^2)*2/(kante.Param.dx*kante.Param.Arho) - kante.Param.A*(kante.KR.y.P-kante.y.P[end])*2/kante.Param.dx - lamda(kante.y.mR,kante.Param.D,kante.Param.A,kante.Param.mu,kante.Param.K)/(2*kante.Param.D*kante.Param.Arho)*abs(kante.y.mR)*kante.y.mR - kante.Param.g*kante.Param.Arho*sin(kante.Param.phi); #-- mR
    dy[k+3*kante.Param.nx+3] = kante.y.eR -(0.5*kante.Param.cv_H2O*(abs(kante.y.mR)*(kante.y.T[end]-kante.KR.y.T)+kante.y.mR*(kante.y.T[end]+kante.KR.y.T)) + kante.Param.A/kante.Param.dx*2*kante.Param.lamW*(kante.y.T[end]-kante.KR.y.T)) #-- eR
end

function lamda(m,D,A,mu,K)
    Re = abs(m)*D/(A*mu); Re = max(Re,1.0e-6)
    lam = 0.25/(log10(K/(3.7*D)+5.74/exp(0.9*log(Re)))^2)
    return lam
end