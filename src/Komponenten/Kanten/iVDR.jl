Base.@kwdef mutable struct iVDR_Param
    I0 = 1
    U0 = 1.5    #!!!Dieser Parameter wird auch auf U_max getestet aber in Matlab nicht!!!
    gamma = 13
    I_max = 1.0e6
end

Base.@kwdef mutable struct y_iVDR
    i::Number = 0.0
end

Base.@kwdef mutable struct iVDR_kante <: Strom_Kante
    #-- geänderte Parameter
    Param::iVDR_Param

    #-- Spannungsknoten links und rechts
    KL::Strom_Knoten
    KR::Strom_Knoten    

    #-- Zustandsvariablen
    y = y_iVDR()

    #-- M-Matrix
    M::Array{Int} = [0]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::iVDR_kante,t)
    #-- Parameter
    (; I0,U0,gamma,I_max) = kante.Param
    #--

    i = kante.y.i

    #-- Spannungsknoten links und rechts
    (; KL,KR,Z) = kante
    UL = KL.y.U
    UR = KR.y.U
    #--

    io = 1.0; if get(Z,"Schaltzeit",0)!=0 io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end

    I = I0*((UL-UR)/U0)^gamma; 
    I = min(max(-I_max,I),I_max);
    dy[k] = i - io*I
end