Base.@kwdef struct WPSP2_Param
    P0 = 1.0e5
    T0 = 293.15
    D = 1.0
    A = pi*(D/2)^2
    cv_H2O = 4182; #-- noch faglich
end

Base.@kwdef mutable struct y_WPSP2
    Param::WPSP2_Param
    M::Number = Param.P0*Param.A/9.81
    T::Number = Param.T0
    P::Number = Param.P0
end

Base.@kwdef mutable struct WPSP2_Knoten <: Wasser_Knoten
    #-- default Parameter
    Param::WPSP2_Param

    #-- Zustandsvariablen 
    y = y_WPSP2(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1; 1; 0] 

    #-- zusätzeliche Infos
    Z::Dict

    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0

    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::WPSP2_Knoten,t)
    #-- Parameter
    (; A,cv_H2O) = knoten.Param
    #--

    (; M, P) = knoten.y

    sum_E = sum_e(knoten.in,knoten.out)
    if typeof(sum_E) == Symbolics.Num
        sum_E = 0.0
    end

    dy[k] =  sum_m(knoten.in,knoten.out)
    dy[k+1] = sum_E/(1e-16*M*cv_H2O)
    dy[k+2] = P-M*9.81/A
end
