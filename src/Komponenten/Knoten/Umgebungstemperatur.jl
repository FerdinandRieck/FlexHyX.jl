Base.@kwdef mutable struct T0_Param
    T0 = 293.15
end

Base.@kwdef mutable struct y_T0
    Param::T0_Param
    T::Number = Param.T0
end

Base.@kwdef mutable struct T0_Knoten <: Knoten
    #-- default Parameter
    Param::T0_Param

    #-- Zustandsvariable
    y = y_T0(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [0] 

    #-- zusätzeliche Infos
    Z::Dict

    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::T0_Knoten,t)
    (; T0) = knoten.Param
    T = knoten.y.T

    (; Z) = knoten

    if haskey(Z,"T0_fcn") 
        dy[k] = T-Z["T0_fcn"](t)
    else
        dy[k] = T-T0
    end 
end

function fcn_T_aussen(t)
    A = 10.0
    P = -A+A*sin(t*2*pi/(24*3600)-pi/2)+273.15
    P = P
end
