Base.@kwdef mutable struct iSP0_Param
    U_soll = 3.6078
    eta = 1.0
end

Base.@kwdef mutable struct y_iSP0
    i_in::Number = 0.0
    i_out::Number = 0.0
end

Base.@kwdef mutable struct iSP0_kante <: Strom_Kante
    #-- geänderte Parameter
    Param::iSP0_Param

    #-- Zustandsvariablen
    y = y_iSP0()

    #-- Spannungsknoten links und rechts
    KL::Strom_Knoten
    KR::Strom_Knoten

    #-- M-Matrix
    M::Array{Int} = [0; 0]

    #-- zusätzliche Infos    
    Z::Dict
end

function Kante!(dy,k,kante::iSP0_kante,t)
    #-- Parameter
    (; U_soll,eta) = kante.Param
    #--

    (; i_in,i_out) = kante.y


    #-- Spannungsknoten links und rechts
    (; KL,KR,Z) = kante
    UL = KL.y.U
    UR = KR.y.U
    #--

    io = 1.0; if get(Z,"Schaltzeit",0)!=0 io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end
    
    dy[k] = ifxaorb(i_in,1/eta,eta)*i_out*UR/UL - i_in;
    dy[k+1] = io*(UR - U_soll) + (1-io)*i_out
end