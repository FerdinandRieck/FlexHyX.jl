Base.@kwdef mutable struct WPHR_Param
    PW0 = 1.0
    T0 = 293.15
    D = 1.0
    A = pi*(D/2)^2
    cv_H2O = 4182; #-- noch faglich
    Ki = 0.01
    Neigung = 1.3
    RT_Soll = 293.15
    Niveau = 0
    T_aussen = 273.15
    ϕ0 = 0.0
    e_max = 0.05
    Jac_init = true
end

Base.@kwdef mutable struct y_WPHR
    Param::WPHR_Param
    M::Number = Param.PW0*1e5*Param.A/9.81
    MT::Number = M*Param.T0 
    P::Number = Param.PW0
    T::Number = Param.T0
    ϕ::Number = Param.ϕ0
    e_zu::Number = 0.0
    E::Number = 0.0
end



Base.@kwdef mutable struct WPHR_Knoten <: Wasser_Knoten
    #-- default Parameter
    Param::WPHR_Param

    #-- Zustandsvariablen 
    y = y_WPHR(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1; 1; 0; 0; 1; 0; 1] 

    #-- zusätzeliche Infos
    Z::Dict

    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0

    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::WPHR_Knoten,t)
    #-- Parameter
    (; A,cv_H2O,Ki,Neigung,RT_Soll,Niveau,e_max,Jac_init) = knoten.Param
    #--

    (; M, MT, P, T, e_zu, ϕ) = knoten.y
    Z = knoten.Z

    io = 1.0
    if Jac_init == true
        io = 1.0
    elseif (haskey(Z,"Schaltzeit")==true) 
        io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) 
    end

    ϕ_max = 1.0; ϕ_min = 0.0; 
    ϕ = min(max(ϕ,ϕ_min),ϕ_max)

    if haskey(Z,"T_aussen") 
        T_aussen = Z["T_aussen"](t)
        DAR = T_aussen - RT_Soll
        VT_Soll = RT_Soll + Niveau - Neigung*DAR*(1.4347 + 0.021*DAR + 247.9*10^-6*DAR^2)
    elseif haskey(Z,"T_soll")
        VT_Soll = Z["T_soll"]
    end

    dy[k] = knoten.sum_m
    dy[k+1] = (knoten.sum_e + e_zu)/(1e-6*cv_H2O)
    dy[k+2] = P-M*9.81/A*1e-5
    dy[k+3] = T-MT/M
    dy[k+4] = Ki*(VT_Soll-T)*ifxaorb((VT_Soll-T),ϕ_max-ϕ,ϕ-ϕ_min)- (1-io)*ϕ
    dy[k+5] = e_zu  - ϕ*e_max
    dy[k+6] = e_zu 
end