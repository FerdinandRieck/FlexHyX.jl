Base.@kwdef mutable struct iS_Param
    R = 1.0
    Jac_init = true
end

Base.@kwdef mutable struct y_iS
    i::Number = 0.0
end

Base.@kwdef mutable struct iS_kante <: Strom_Kante
    #-- geänderte Parameter
    Param::iS_Param

    #-- Spannungsknoten links und rechts
    KL::Strom_Knoten
    KR::Strom_Knoten    

    #-- Zustandsvariablen
    y = y_iS()

    #-- M-Matrix
    M::Array{Int} = [0]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::iS_kante,t)
    #-- Parameter
    (; R,Jac_init) = kante.Param
    #--

    i = kante.y.i

    #-- Spannungsknoten links und rechts
    (; KL,KR,Z) = kante
    UL = KL.y.U
    UR = KR.y.U
    #--

    io = 1.0
    if Jac_init == true
        io = 1.0
    elseif (haskey(Z,"Schaltzeit")==true) 
        io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) 
    end
   
    dy[k] = i - io*(UL-UR)/R;
end