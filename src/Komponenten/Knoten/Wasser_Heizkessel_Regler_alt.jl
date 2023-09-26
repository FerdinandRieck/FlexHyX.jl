Base.@kwdef struct WPHR_Param
    P0 = 1.0e5
    T0 = 293.15
    D = 1.0
    A = pi*(D/2)^2
    cv_H2O = 4182; #-- noch faglich
    Kp = 1000
    Neigung = 1.3
    RT_Soll = 293.15
    Niveau = 1
    T_aussen = 273.15
end

Base.@kwdef mutable struct y_WPHR
    Param::WPHR_Param
    M::Number = Param.P0*Param.A/9.81
    MT::Number = M*Param.T0 
    P::Number = Param.P0
    T::Number = Param.T0
    e_zu::Number = 0.0
    VT_Soll::Number = 0.0
end



Base.@kwdef mutable struct WPHR_Knoten <: Wasser_Knoten
    #-- default Parameter
    Param::WPHR_Param

    #-- Zustandsvariablen 
    y = y_WPHR(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1; 1; 0; 0; 1; 0] 

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
    (; A,cv_H2O,Kp,Neigung,RT_Soll,Niveau) = knoten.Param
    #--

    (; M, MT, P, T, e_zu, VT_Soll) = knoten.y
    Z = knoten.Z



    dy[k] =  sum_m(knoten.in,knoten.out)
    dy[k+1] = (knoten.sum_e + e_zu)/cv_H2O
    dy[k+2] = P-M*9.81/A
    dy[k+3] = T-MT/M
    dy[k+4] = Kp*(VT_Soll-T) 
    if haskey(Z,"T_aussen") 
        T_aussen = Z["T_aussen"](t)
        DAR = T_aussen - RT_Soll
        dy[k+5] = VT_Soll - (RT_Soll * Niveau - Neigung*DAR*(1.4347 + 0.021*DAR + 247.9*10^-6*DAR^2))
    elseif haskey(Z,"T_soll")
        dy[k+5] = VT_Soll - Z["T_soll"]
    end
end