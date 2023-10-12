Base.@kwdef mutable struct WP2_Param
    cv_H2O = 4182.0
    P0 = 1.0e5
    T0 = 293.15
end

Base.@kwdef mutable struct y_WP2
    Param::WP2_Param
    P::Number = Param.P0
    T::Number = Param.T0
end

Base.@kwdef mutable struct WP2_Knoten <: Wasser_Knoten
    #-- geänderte Parameter
    Param::WP2_Param

    #-- Zustandsvariablen
    y = y_WP2(Param=Param)
    
    #-- M-Matrix
    M::Array{Int} = [0; 1] 

    #-- zusätzeliche Infos
    Z::Dict
    
    #-- Knotenbilanz
    sum_m::Number = 0.0
    sum_e::Number = 0.0   
    
    #-- Kanten
    in::Array{Any} = []
    out::Array{Any} = []
end

function Knoten!(dy,k,knoten::WP2_Knoten,t)
    (; cv_H2O) = knoten.Param

    sum_E = sum_e(knoten.in,knoten.out)
    if typeof(sum_E) == Symbolics.Num
        sum_E = 0.0
    end

    #sum_E = -abs(max(min(1e-6,sum_E),-1e-6))
    if t > 980
        if knoten.Z["Nr"] == 15
           # @show sum_E
        end
    end

   

    Masse = 0.1
    dy[k] = sum_m(knoten.in,knoten.out)
    dy[k+1] = sum_E/(1e-6*Masse*cv_H2O) 
end

#=
function Knoten!(dy,k,knoten::WP2_Knoten,t)
    (; cv_H2O) = knoten.Param

    Masse = 0.1
    dy[k] = sum_m(knoten.in,knoten.out)
    dy[k+1] =  knoten.sum_e/(Masse*cv_H2O*1e-16) 
end
=#