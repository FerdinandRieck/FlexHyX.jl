Base.@kwdef mutable struct iDL_Param
    R = 0.01
    alpha = 1
    L0 = 1/R
end

Base.@kwdef mutable struct y_iDL
    Param::iDL_Param
    i::Number = 0.0
    L::Number = Param.L0
end

Base.@kwdef mutable struct iDL_kante <: Strom_Kante
    #-- geänderte Parameter
    Param::iDL_Param

    #-- Zustandsvariablen
    y = y_iDL(Param=Param)

    #-- Spannungsknoten links und rechts
    KL::Strom_Knoten
    KR::Strom_Knoten

    M::Array{Int} = [0; 1]
    Z::Dict
end

function Kante!(dy,k,kante::iDL_kante,t)
    #-- Parameter
    (; R,alpha) = kante.Param
    #--

    iDL = kante.y.i
    L = kante.y.L

    #-- Spannungsknoten links und rechts
    (; KL,KR,Z) = kante
    UL = KL.y.U
    UR = KR.y.U
    #--

    io = 1.0; if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end
    U = UL-UR
    L_soll = 1/R
    L1=min(max(L,0),L_soll)
    res = U/alpha*ifxaorb(U,L_soll-L1,L1)
    dy[k] = iDL - io*(UL-UR)*L
    dy[k+1] = io*res
end