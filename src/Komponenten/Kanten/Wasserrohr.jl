Base.@kwdef mutable struct mWRo_Param
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
    WENO = true
    fluxPL = Array{Number}(undef, nx+1)
    fluxPR = Array{Number}(undef, nx+1)
    fluxmL = Array{Number}(undef, nx+1)
    fluxmR = Array{Number}(undef, nx+1)
    fluxTL = Array{Number}(undef, nx+1)
    fluxTR = Array{Number}(undef, nx+1)
end

Base.@kwdef mutable struct y_mWRo
    mL::Number = 0.0
    eL::Number = 0.0
    P
    _m
    T
    mR::Number = 0.0
    eR::Number = 0.0
end

Base.@kwdef mutable struct mWRo_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWRo_Param

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- Zustandsvariablen
    y = y_mWRo(P = f(Vector(1:2:2*Param.nx), KL.y.P, KR.y.P, Param.nx), 
               T = f(Vector(1:2:2*Param.nx), KL.y.T, KR.y.T, Param.nx),
               _m = zeros(Param.nx)
               )

    #-- M-Matrix
    M::Array{Int} = [1; 0; ones(Int,3*Param.nx); 1; 0]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWRo_kante,t)
    #-- Parameter
    (; nx,dx,a2,leit,Arho,A,D,cv_H2O,mu,K,lamW,phi,g,WENO,fluxPL,fluxPR,fluxmL,fluxmR,fluxTL,fluxTR) = kante.Param
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
    dy[k] = -(m[1]^2-mL^2)*2/(dx*Arho) - A*(P[1]-PL)*2/dx - lambda(mL,D,A,mu,K)/(2*D*Arho)*abs(mL)*mL - g*Arho*sin(phi); #-- mL
    dy[k+1] = eL -(0.5*cv_H2O*(abs(mL)*(TL-T[1])+mL*(TL+T[1])) + A/dx*2*lamW*(TL-T[1])); #-- eL
    
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
        dy[k+i+1] = -a2/A*(fRm-fLm)/dx
        dy[k+i+1+nx] = -(fRm^2-fLm^2)/(dx*Arho) - A*(fRP-fLP)/dx - lambda(m[i],D,A,mu,K)/(2*D*Arho)*abs(m[i])*m[i] - g*Arho*sin(phi)
        dy[k+i+1+nx*2] = -1/Arho*m[i]*ifxaorb(m[i],T[i]-fLT,fRT-T[i])*2/dx + leit*2/(dx^2)*(fLT-2*T[i]+fRT)
    end
    
    #-- Rohr rechts
    dy[k+3*nx+2] = -(mR^2-m[end]^2)*2/(dx*Arho) - A*(PR-P[end])*2/dx - lambda(mR,D,A,mu,K)/(2*D*Arho)*abs(mR)*mR - g*Arho*sin(phi); #-- mR
    dy[k+3*nx+3] = eR -(0.5*cv_H2O*(abs(mR)*(T[end]-TR)+mR*(T[end]+TR)) + A/dx*2*lamW*(T[end]-TR)) #-- eR
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


#= Benchmark Tests mit revover_weno Funktionen
#--- out-of-place nx = 10
BenchmarkTools.Trial: 4 samples with 1 evaluation.
 Range (min … max):  1.371 s …    1.638 s  ┊ GC (min … max): 4.10% … 7.75%
 Time  (median):     1.477 s               ┊ GC (median):    3.63%
 Time  (mean ± σ):   1.491 s ± 123.299 ms  ┊ GC (mean ± σ):  4.76% ± 2.09%

  █       █                            █                   █  
  █▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.37 s         Histogram: frequency by time         1.64 s <

 Memory estimate: 320.46 MiB, allocs estimate: 6799535.

 #--- in-place nx = 10
 BenchmarkTools.Trial: 3 samples with 1 evaluation.
 Range (min … max):  1.466 s …    2.092 s  ┊ GC (min … max): 3.38% … 3.56%
 Time  (median):     1.606 s               ┊ GC (median):    3.45%
 Time  (mean ± σ):   1.721 s ± 328.624 ms  ┊ GC (mean ± σ):  3.48% ± 0.09%

  █            █                                           █  
  █▁▁▁▁▁▁▁▁▁▁▁▁█▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  1.47 s         Histogram: frequency by time         2.09 s <

 Memory estimate: 296.10 MiB, allocs estimate: 6623482.

 #--- in-place nx = 100
 BenchmarkTools.Trial: 1 sample with 1 evaluation.
 Single result which took 39.846 s (12.59% GC) to evaluate,
 with a memory estimate of 23.43 GiB, over 343769246 allocations.

 #--- out-of-place nx = 100
 BenchmarkTools.Trial: 1 sample with 1 evaluation.
 Single result which took 47.508 s (13.88% GC) to evaluate,
 with a memory estimate of 27.67 GiB, over 359363538 allocations.
=#