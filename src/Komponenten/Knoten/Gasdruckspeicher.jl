Base.@kwdef mutable struct GPSP_Param
    P0 = 144210
    T0 = 293.15
    V = 1.0
    Rs = 4124.2
    cv = 1.01798*1.0e4
end

Base.@kwdef mutable struct y_GPSP
    Param::GPSP_Param
    M::Number = Param.P0*Param.V/(Param.T0*Param.Rs)
    MT::Number = M*Param.T0 
    P::Number = Param.P0
    T::Number = Param.T0
end

Base.@kwdef mutable struct GPSP_Knoten <: Gas_Knoten
    #-- default Parameter
    Param::GPSP_Param

    #-- Zustandsvariablen 
    y = y_GPSP(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1; 1; 0; 0] 

    #-- Jacobi Struktur
    J::Array{Int} = [0 0 0 0; 0 0 0 0; 0 1 1 0; 1 1 0 1]
    J_fluss::Array{String} = ["sum_m"; "sum_e"]

    #-- zusätzeliche Infos
    Z::Dict
end

function Knoten!(dy,k,sum_i,sum_m,sum_e,knoten::GPSP_Knoten)
    #-- Parameter
    (; V,Rs,cv) = knoten.Param
    #--
    (; M, MT, P, T) = knoten.y
    dy[k] = sum_m
    dy[k+1] = sum_e/cv
    dy[k+2] = P-Rs*MT/V
    dy[k+3] = T-MT/max(M,1.0e-9)
end