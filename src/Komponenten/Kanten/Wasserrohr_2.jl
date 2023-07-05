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
    P_vec = PL:(PR-PL)*0.5/nx:PR
    T_vec = TL:(TR-TL)*0.5/nx:TR
    P = Vector(P_vec[2:2:end-1])
    T = Vector(T_vec[2:2:end-1])
end

#-- Rohrränder ---------------------------------
Base.@kwdef mutable struct y_mWRo2
    Param::mWRo2_Param
    mL::Number = 0.0
    eL::Number = 0.0
    P::Array{Number} = Param.P
    _m::Array{Number} = zeros(Param.nx)
    T::Array{Number} = Param.T
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
    (; nx,dx,a2,leit,Arho,A,D,cv_H2O,mu,K,lamW,phi,g) = kante.Param
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

    (; KL,KR) = kante
    PL = KL.y.P
    TL = KL.y.T
    PR = KR.y.P
    TR = KR.y.T

   # P = [PL; y[k+2:k+1+nx]; PR]
   # m = [mL; y[k+2+nx:k+1+2*nx]; mR]
   # T = [TL; y[k+2+2*nx:k+1+3*nx]; TR]

   # Re = abs(m)*D/(A*var.mu); Re = max(Re,1.0e-6);
   # lam = 0.25./(log10(K/(3.7*D)+5.74./exp(0.9*log(Re))).^2);

    #-- Rohr links
    dy[k] = -(m[1]^2-mL^2)*2/(dx*Arho) - A*(P[1]-PL)*2/dx - lamda(mL,D,A,mu,K)/(2*D*Arho)*abs(mL)*mL - g*Arho*sin(phi); #-- mL
    dy[k+1] = eL -(0.5*cv_H2O*(abs(mL)*(TL-T[1])+mL*(TL+T[1])) + A/dx*2*lamW*(TL-T[1])); #-- eL
    #-- Rohr mitte
    for i = 1:nx
        if i==1
            P_m12 = PL; m_m12 = mL; T_m12 = TL;
        else
            P_m12 = 0.5*(P[i]+P[i-1]); m_m12 = 0.5*(m[i]+m[i-1]); T_m12 = 0.5*(T[i]+T[i-1]);
        end
        @show P_m12
        if i==nx
            P_p12 = PR; m_p12 = mR; T_p12 = TR;
        else
            P_p12 = 0.5*(P[i]+P[i+1]); m_p12 = 0.5*(m[i]+m[i+1]); T_p12 = 0.5*(T[i]+T[i+1]);
        end
        dy[k+i*3-1] = -a2/A*(m_p12-m_m12)/dx
        dy[k+i*3] = -(m_p12^2-m_m12^2)/(dx*Arho) - A*(P_p12-P_m12)/dx - lamda(m[i],D,A,mu,K)/(2*D*Arho)*abs(m[i])*m[i] - g*Arho*sin(phi)
        dy[k+i*3+1] = -1/Arho*m[i]*ifxaorb(m[i],T[i]-T_m12,T_p12-T[i])*2/dx + leit*2/(dx^2)*(T_m12-2*T[i]+T_p12)
    end
    #-- Rohr rechts
    dy[k+3*nx+1] = -(mR^2-m[end]^2)*2/(dx*Arho) - A*(PR-P[end])*2/dx - lamda(mR,D,A,mu,K)/(2*D*Arho)*abs(mR)*mR - g*Arho*sin(phi); #-- mR
    dy[k+3*nx+2] = eR -(0.5*cv_H2O*(abs(mR)*(T[end]-TR)+mR*(T[end]+TR)) + A/dx*2*lamW*(T[end]-TR)) #-- eR
end

function lamda(m,D,A,mu,K)
    Re = abs(m)*D/(A*mu); Re = max(Re,1.0e-6)
    lam = 0.25/(log10(K/(3.7*D)+5.74/exp(0.9*log(Re)))^2)
    return lam
end