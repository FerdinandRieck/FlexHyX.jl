Base.@kwdef mutable struct mWTaRM_Param
    nx = 1
    L = 0.286
    dx = L/max(nx,1/2)
    rho0 = 1000
    Di = 0.02
    Dm = 0.022
    Da = 0.1
    Aa = pi/4*(Da^2-Dm^2)
    Ai = pi/4*Di^2
    Aarho = Aa*rho0
    Airho = Ai*rho0
    lamW = 0.6
    lamRohr = 401.0 #-- Wärmleitung Kupfer
    kA = 2000.0
    cv_H2O = 4182.0; #-- noch faglich
    leit = lamW/(rho0*cv_H2O)
    mu = 1.0e-3
    WENO = true
    Richtung = "gleich"
    Ringspalt = true
    fluxTL = Array{Number}(undef, nx+1)
    fluxTR = Array{Number}(undef, nx+1)
    m_dot = 1.0
end

Base.@kwdef mutable struct y_mWTaRM
    m::Number = 0.0
    eL::Number = 0.0
    eR::Number = 0.0
    T
end


Base.@kwdef mutable struct mWTaRM_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWTaRM_Param

    #-- Wasserknoten links und rechts
    KL::Wasser_Knoten
    KR::Wasser_Knoten

    #-- Zustandsvariablen
    y = y_mWTaRM(T = f(Vector(1:2:2*Param.nx), KL.y.T, KR.y.T, Param.nx))

    #-- Ringspalt/Rohr innen
    R = 0

    #-- M-Matrix
    M::Array{Int} = [0; 0; 0; ones(Int,Param.nx)]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWTaRM_kante,t)
    #-- Parameter
    (; nx,dx,leit,Ai,Aa,Aarho,Airho,Dm,cv_H2O,lamW,WENO,Richtung,fluxTL,fluxTR,Ringspalt) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    eL = kante.y.eL
    eR = kante.y.eR
    T = kante.y.T
    #--

    (; KL,KR,R,Z) = kante
    TL = KL.y.T
    TR = KR.y.T
    mL = KL.in[1].y.m
    T_aussen = R.y.T

    if WENO == true
        recover_weno!(T,fluxTL,fluxTR) #??? hier vieleicht nur recover!(T), da bereits ubwind Diskertisierung mit ifxaorb(m[i],T[i]-fLT,fRT-T[i])??? 
    else
        recover!(T,fluxTL,fluxTR)
    end

    if (Richtung == "gegen") && nx > 1
        T_aussen = reverse(T_aussen)
    end

    if Ringspalt==true
        if isa(Z["kA"],Function)
            kA = Z["kA"](kante)
        else
            kA = Z["kA"]
        end
        Arho = Aarho
        A = Aa
    else
        kA = Z["kA"]
        Arho = Airho
        A = Ai
    end

    #-- Rohr links
    dy[k] = m - mL #-- m
    TRL = T[1] - (T[1]-T[2])/dx * -0.5*dx
    dy[k+1] = eL -1e-6*(0.5*cv_H2O*(abs(m)*(TL-TRL)+m*(TL+TRL)) + A/dx*2*lamW*(TL-TRL)) #-- eL
    #-- Rohr rechts
    TRR = T[nx-1] - (T[nx-1]-T[nx])/dx * 1.5*dx
    dy[k+2] = eR -1e-6*(0.5*cv_H2O*(abs(m)*(TRR-TR)+m*(TRR+TR)) + A/dx*2*lamW*(TRR-TR)) #-- eR

    #-- Rohr mitte
    fRT = TL #-- linke Randbedingung
    for i = 1:nx
        fLT = fRT
        if i == nx
            fRT = TR
        else 
            fRT = 0.5*(fluxTL[i+1]+fluxTR[i+1])
        end
        dy[k+i+2] = -1/Arho*m*ifxaorb(m,T[i]-fLT,fRT-T[i])*2/dx + leit*2/(dx^2)*(fLT-2*T[i]+fRT) - kA*pi*Dm/(Arho*cv_H2O)*(T[i]-T_aussen[i])  #-- T
    end
end

function mWTaRM_init(knoten,kanten,M,kk,von,nach)
    Params = MakeParam(kk) 
    kante = mWTaRM_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk)
    push!(kanten,kante)
    append!(M, kante.M)
    if haskey(kk,"RefKante") 
        if isa(kk["RefKante"],Int)
            R = kk["RefKante"]
            kanten[end].R = kanten[R]
            kanten[R].R = kanten[end]
        end
    end
end