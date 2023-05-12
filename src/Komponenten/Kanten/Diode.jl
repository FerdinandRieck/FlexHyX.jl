Base.@kwdef mutable struct iD_Param
    R = 0.01
    I0 = 0.001
    U0 = 0.5
    f = beta -> [I0*beta[1]*exp(beta[1]*U0) - 1/R]
    res = nlsolve(f,[0.1]); 
    beta = res.zero[1]
    a = -U0/R+I0*(exp(beta*U0)-1);
end

Base.@kwdef mutable struct y_iD
    i::Number = 0.0
end

Base.@kwdef mutable struct iD_kante <: Strom_Kante
    #-- geänderte Parameter
    Param::iD_Param

    #-- Zustandsvariablen
    y = y_iD()

    #-- Spannungsknoten links und rechts
    KL::Strom_Knoten
    KR::Strom_Knoten

    M::Array{Int} = [0]
    Z::Dict
end

function Kante!(dy,k,kante::iD_kante,t)
    #-- Parameter
    (; R,I0,U0,beta) = kante.Param
    #--

    iD = kante.y.i

    #-- Spannungsknoten links und rechts
    (; KL,KR,Z) = kante
    UL = KL.y.U
    UR = KR.y.U
    #--

    io = 1.0; if get(Z,"Schaltzeit",0)!=0 io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end
    U = UL-UR
    I = I0*(exp(beta*min(U,U0))-1) + max(U-U0,0)/R;
    
    dy[k] = iD - io*I;
end