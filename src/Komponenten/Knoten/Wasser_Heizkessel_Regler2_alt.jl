Base.@kwdef struct WPHR2_Param
    P0 = 1.0e5
    T0 = 293.15
    D = 1.0
    A = pi*(D/2)^2
    cv_H2O = 4182; #-- noch faglich
    Kp = 0.01
    Neigung = 1.3
    RT_Soll = 293.15
    Niveau = 1
    T_aussen = 273.15
    ϕ0 = 0.5
    e_max = 5.0e4
end

Base.@kwdef mutable struct y_WPHR2
    Param::WPHR2_Param
    M::Number = Param.P0*Param.A/9.81
    MT::Number = M*Param.T0 
    P::Number = Param.P0
    T::Number = Param.T0
    ϕ::Number = Param.ϕ0
    e_zu::Number = 0.0
    E::Number = 0.0
    VT_Soll::Number = 0.0
end



Base.@kwdef mutable struct WPHR2_Knoten <: Wasser_Knoten
    #-- default Parameter
    Param::WPHR2_Param

    #-- Zustandsvariablen 
    y = y_WPHR2(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1; 1; 0; 0; 1; 0; 1; 0] 

    #-- zusätzeliche Infos
    Z::Dict

    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0

    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::WPHR2_Knoten,t)
    #-- Parameter
    (; A,cv_H2O,Kp,Neigung,RT_Soll,Niveau,e_max) = knoten.Param
    #--

    (; M, MT, P, T, e_zu, ϕ, VT_Soll) = knoten.y
    Z = knoten.Z

    io = 1.0;  
    if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end

    ϕ_max = 1.0; ϕ_min = 0.0; 
    ϕ = min(max(ϕ,ϕ_min),ϕ_max)

    dy[k] = sum_m(knoten.in,knoten.out)
    dy[k+1] = (knoten.sum_e + e_zu)/cv_H2O
    dy[k+2] = P-M*9.81/A
    dy[k+3] = T-MT/M
    dy[k+4] = Kp*(VT_Soll-T)*ifxaorb((VT_Soll-T),ϕ_max-ϕ,ϕ-ϕ_min) - (1-io)*ϕ
    dy[k+5] = e_zu  - ϕ*e_max
    dy[k+6] = e_zu 
    if haskey(Z,"T_aussen") 
        T_aussen = Z["T_aussen"](t)
        DAR = T_aussen - RT_Soll
        dy[k+7] = VT_Soll - (RT_Soll * Niveau - Neigung*DAR*(1.4347 + 0.021*DAR + 247.9*10^-6*DAR^2))
    elseif haskey(Z,"T_soll")
        dy[k+7] = VT_Soll - Z["T_soll"]
    end
end