Base.@kwdef mutable struct eWTL_Param
   cL = 1006.7
   t_wechsel = 3*3600
   t_end = 700000
   t_lüften = 300
   Jac_init = true
end

Base.@kwdef mutable struct y_eWTL
    e::Number = 0.0
end

Base.@kwdef mutable struct eWTL_kante <: Temp_Kante
    #-- default Parameter
    Param::eWTL_Param

    #-- Knoten links und rechts
    KL::Knoten  
    KR::Knoten

    #-- Zustandsvariablen
    y = y_eWTL()

    #-- M-Matrix
    M::Array{Int} = [0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::eWTL_kante,t)
    #-- Parameter
    (; cL,t_wechsel,t_end,Jac_init,t_lüften) = kante.Param
    #--

    #-- Zustandsvariable
    e = kante.y.e
    #--

    (; KL,KR,Z) = kante
    TL = KL.y.T
    TR = KR.y.T

    io = 1.0
    if Jac_init == true
        io = 1.0
    elseif (haskey(Z,"Schaltzeit")==true) 
        io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) 
        t_lüften = Z["Schaltzeit"][end] - Z["Schaltzeit"][end-1]
    else #-- beim ersten Druchgang Schaltzeiten für Luftwechsel berechnen wenn nicht angegeben
        n = t_end/t_wechsel
        n = trunc(Int, n)*2+1   #-- Anzahl Schaltzeiten muss ungerade sein
        Schaltzeit = zeros(n,1)
        Schaltzeit[1] = 50000
        for i = 2:2:n-1
            Schaltzeit[i] = Schaltzeit[i-1] + t_lüften
            Schaltzeit[i+1] = Schaltzeit[i] + t_wechsel 
        end
        Z["Schaltzeit"] = Schaltzeit
        io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"])
    end

    dot_mL = KL.Param.Masse/t_lüften

    dy[k] = e - 1e-6*io*dot_mL*cL*(TL-TR)
end