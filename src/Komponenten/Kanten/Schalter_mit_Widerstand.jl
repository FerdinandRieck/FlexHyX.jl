Base.@kwdef mutable struct iS_Param
    R = 1.0
end

Base.@kwdef mutable struct y_iS
    i::Number = 0.0
end

Base.@kwdef mutable struct iS_kante <: Strom_Kante
    #-- geänderte Parameter
    Param::iS_Param

    #-- Zustandsvariablen
    y = y_iS()

    #-- Spannungsknoten links und rechts
    KL::Strom_Knoten
    KR::Strom_Knoten

    #-- M-Matrix
    M::Array{Int} = [0]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::iS_kante,t)
    #-- Parameter
    (; R) = kante.Param
    #--

    i = kante.y.i

    #-- Spannungsknoten links und rechts
    (; KL,KR,Z) = kante
    UL = KL.y.U
    UR = KR.y.U
    #--

    io = 1.0; if get(Z,"Schaltzeit",0)!=0 io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end
    dy[k] = i - io*(UL-UR)/R;
end