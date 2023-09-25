Base.@kwdef mutable struct eWTL_Param
   cL = 1006.7
end

Base.@kwdef mutable struct y_eWTL
    e::Number = 0.0
end

Base.@kwdef mutable struct eWTL_kante <: Temp_Kante
    #-- default Parameter
    Param::eWTL_Param

    #-- Zustandsvariablen
    y = y_eWTL()

    #-- Gasknoten links und rechts
    KL::Knoten  
    KR::Knoten

    #-- M-Matrix
    M::Array{Int} = [0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::eWTL_kante,t)
    #-- Parameter
    (; cL) = kante.Param
    #--

    #-- Zustandsvariable
    e = kante.y.e
    #--

    (; KL,KR,Z) = kante
    TL = KL.y.T
    TR = KR.y.T

    io = 1.0
    if (haskey(Z,"Schaltzeit")==true) 
        io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) 
        t_lüften = Z["Schaltzeit"][end] - Z["Schaltzeit"][end-1]
    else
        t_lüften = 3600
    end

    M = KL.Param.Masse/t_lüften

    dy[k] = e - io*M*cL*(TL-TR) 
end