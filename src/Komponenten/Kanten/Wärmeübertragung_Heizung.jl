Base.@kwdef mutable struct mWTH_Param
    nx = 1
    L = 1.0
    dx = L/max(nx,1/2)
    rho0 = 1000.0
    M = 1.0
    A = M/(rho0*L)
    Aü = 1.89
    lam = 0.6
    cv = 4182.0; #-- noch faglich
    Arho = A*rho0
    leit = lam/(rho0*cv)
    kA = 380.0
    WENO = true
    fluxTL = Array{Number}(undef, nx+1)
    fluxTR = Array{Number}(undef, nx+1)
    m_dot = 1.0
end

Base.@kwdef mutable struct y_mWTH
    m::Number = 0.0
    eL::Number = 0.0
    eR::Number = 0.0
    T_Raum::Number = 0.0
    T_Aussen::Number = 0.0
    T
end


Base.@kwdef mutable struct mWTH_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWTH_Param

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- Zustandsvariablen
    y = y_mWTH(T = f(Vector(1:2:2*Param.nx), KL.y.T, KR.y.T, Param.nx))

    #-- M-Matrix
    M::Array{Int} = [0; 0; 0; 1; 1; ones(Int,Param.nx)]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWTH_kante,t)
    #-- Parameter
    (; nx,dx,leit,Arho,A,Aü,M,cv,lam,kA,WENO,fluxTL,fluxTR) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    eL = kante.y.eL
    eR = kante.y.eR
    T_Raum = kante.y.T_Raum
    T_Aussen = kante.y.T_Aussen
    T = kante.y.T
    #--

    (; KL,KR,Z) = kante
    TL = KL.y.T
    TR = KR.y.T
    mL = KL.in[1].y.m

    if WENO == true
        recover_weno!(T,fluxTL,fluxTR) #??? hier vieleicht nur recover!(T), da bereits ubwind Diskertisierung mit ifxaorb(m[i],T[i]-fLT,fRT-T[i])??? 
    else
        recover!(T,fluxTL,fluxTR)
    end

    io = 1.0; if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end    

    #-- Rohr links
    dy[k] = m - mL #-- m
    dy[k+1] = eL - m*(TL-T[1]) #- A/dx*2*lamW*(TL-T[1]) Wärmeübetragung weglassen
    #dy[k+1] = eL -(0.5*cv*(abs(m)*(TL-T[1])+m*(TL+T[1])) + A/dx*2*lam*(TL-T[1])); #-- eL
    #-- Rohr rechts
    dy[k+2] = eR - m*(T[nx]-TR) #- A/dx*2*lamW*(T[nx]-TR) Wärmeübertragung weglassen
    #dy[k+2] = eR -(0.5*cv*(abs(m)*(T[end]-TR)+m*(T[end]+TR)) + A/dx*2*lam*(T[end]-TR)) #-- eR
    if haskey(Z,"Q_dot") 
        if isa(Z["Q_dot"],Number) dy[k+2] = dy[k+2] - io*Z["Q_dot"]/cv end 
        if isa(Z["Q_dot"],Function) dy[k+2] = dy[k+2] - io*Z["Q_dot"](t,kante) end
    end 
    dy[k+3] = 
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
        if haskey(Z,"T_aussen")
            dy[k+i+2] = dy[k+i+2] - kA*Aü/(M*cv)*(T[i]-Z["T_aussen"])
        end
        if haskey(Z,"KnotenAussen")
            T_aussen = KA.y.T
            dy[k+i+2] = dy[k+i+2] - kA*Aü/(M*cv)*(T[i]-T_aussen)
        end
    end
end

function fcn_Q_dot2(t,kante)
    A = 1000.0
    P = 2*A+A*sin(t*2*pi/(24*3600))
    P = -P/kante.Param.cv
end

function fcn_Q_dot(t,kante)
    #m = kante.KL.in[1].KL.in[1].y.m
    #TL = kante.KL.in[1].KL.in[1].KL.y.T
    #TR = kante.KL.in[1].KL.in[1].KR.y.T
    #@show  kante.KL.in[1].KL.in[1].KL.y.T
    m = kante.y.m
    TL = kante.KL.y.T
    TR = kante.KR.y.T
    e = m*(TL-TR)
    return e*-1.0
end