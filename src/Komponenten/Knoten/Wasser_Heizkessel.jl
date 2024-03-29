Base.@kwdef mutable struct WPH_Param
    PW0 = 1.0
    T0 = 293.15
    D = 1.0
    A = pi*(D/2)^2
    cv_H2O = 4182; #-- noch faglich
    Jac_init = true
end

Base.@kwdef mutable struct y_WPH
    Param::WPH_Param
    M::Number = Param.PW0*1e5*Param.A/9.81
    MT::Number = M*Param.T0 
    P::Number = Param.PW0
    T::Number = Param.T0
    e_zu::Number = 0.0
end



Base.@kwdef mutable struct WPH_Knoten <: Wasser_Knoten
    #-- default Parameter
    Param::WPH_Param

    #-- Zustandsvariablen 
    y = y_WPH(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1; 1; 0; 0; 0] 

    #-- zusätzeliche Infos
    Z::Dict

    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0

    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::WPH_Knoten,t)
    #-- Parameter
    (; A,cv_H2O,Jac_init) = knoten.Param
    #--

    (; M, MT, P, T, e_zu) = knoten.y
    Z = knoten.Z

    io = 1.0
    if Jac_init == true
        io = 1.0
    elseif (haskey(Z,"Schaltzeit")==true) 
        io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) 
    end

    dy[k] =  sum_m(knoten.in,knoten.out)
    dy[k+1] = (knoten.sum_e + e_zu)/(1e-6*cv_H2O)
    dy[k+2] = P-M*9.81/A*1e-5
    dy[k+3] = T-MT/M
    dy[k+4] = e_zu + io*knoten.sum_e*1.1
end