Base.@kwdef struct WPSP3_Param
    P0 = 1.0e5
    T0 = 293.15
    D = 1.0
    A = pi*(D/2)^2
    cv_H2O = 4182; #-- noch faglich
end

Base.@kwdef mutable struct y_WPSP3
    Param::WPSP3_Param
    M::Number = Param.P0*1e5*Param.A/9.81
    MT::Number = M*Param.T0 
    P::Number = Param.P0
    T::Number = Param.T0
end

Base.@kwdef mutable struct WPSP3_Knoten <: Wasser_Knoten
    #-- default Parameter
    Param::WPSP3_Param

    #-- Zustandsvariablen 
    y = y_WPSP3(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1; 1; 0; 0] 

    #-- zusätzeliche Infos
    Z::Dict

    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0

    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::WPSP3_Knoten,t)
    #-- Parameter
    (; A,cv_H2O) = knoten.Param
    #--

    (; M, MT, P, T) = knoten.y

    sum_E = sum_e(knoten.in,knoten.out)
    if typeof(sum_E) == Symbolics.Num
        sum_E = 0.0
    end



    dy[k] = knoten.sum_m
    dy[k+1] = sum_E/(cv_H2O*1e-6)
    dy[k+2] = P-M*9.81/A*1e-5
    dy[k+3] = T-MT/M
end

#=
function Knoten!(dy,k,knoten::WPSP3_Knoten,t)
    #-- Parameter
    (; A,cv_H2O) = knoten.Param
    #--

    (; M, MT, P, T) = knoten.y

    dy[k] = knoten.sum_m
    dy[k+1] = knoten.sum_e/(1e-16*cv_H2O)
    dy[k+2] = P-M*9.81/A*1e-5
    dy[k+3] = T-MT/M
end
=#