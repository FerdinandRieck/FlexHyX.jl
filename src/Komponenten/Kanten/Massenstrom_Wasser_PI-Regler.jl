Base.@kwdef mutable struct mWfPI_Param
    Kp = 10.0
    Ki = 1.0
    T_soll = 293.15
    cv_H2O = 4182.0
    Jac_init = true
end

Base.@kwdef mutable struct y_mWfPI
    Param::mWfPI_Param
    m::Number = 0.0
    e::Number = 0.0
    T_err_dt::Number = 0.0
end

Base.@kwdef mutable struct mWfPI_kante <: Wasser_Kante
    #-- default Parameter
    Param::mWfPI_Param

    #-- Wasserknoten links und rechts
    KL::Knoten
    KR::Knoten
    K = [] #-- Knoten vom Haus

    #-- Zustandsvariablen
    y = y_mWfPI(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [0; 0; 1] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWfPI_kante,t)
    #-- Parameter
    (; T_soll, Kp, Ki, cv_H2O ,Jac_init) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    e = kante.y.e
    T_err_dt = kante.y.T_err_dt
    #--

    (; KL,KR,K,Z) = kante
    TL = KL.y.T
    TR = KR.y.T
    T = K.y.T #-- Zimmertemperatur im Haus

    io = 1.0
    if Jac_init == true
        io = 1.0
    elseif (haskey(Z,"Schaltzeit")==true) 
        io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) 
    end

    u = Kp*(T_soll-T) + Ki*T_err_dt

    dy[k] =  m - u
    dy[k+1] = e - 1e-6*cv_H2O*m*ifxaorb(m,TL,TR)
    dy[k+2] = T_soll-T
end

function mWfPI_init(knoten,kanten,M,kk,von,nach)
    Params = MakeParam(kk) 
    kante = mWfPI_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk) 
    push!(kanten,kante)
    append!(M, kante.M)
    K = kk["RefKnoten"]
    kanten[end].K = knoten[K]
end