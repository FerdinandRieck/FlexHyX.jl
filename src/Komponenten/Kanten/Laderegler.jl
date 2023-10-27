Base.@kwdef mutable struct iLR_Param
    R = 0.01
    alpha = 10
    P_max = 500
    U_soll = 13.5
    i_max = 10 #???
    Jac_init = true
end

Base.@kwdef mutable struct y_iLR
    i::Number = 0.0
end

Base.@kwdef mutable struct iLR_kante <: Strom_Kante
    #-- geänderte Parameter
    Param::iLR_Param

    #-- Spannungsknoten links und rechts
    KL::Strom_Knoten
    KR::Strom_Knoten

    #-- Zustandsvariablen
    y = y_iLR()

    #-- M-Matrix
    M::Array{Int} = [1]

    #-- zusätzliche Infos    
    Z::Dict
end

function Kante!(dy,k,kante::iLR_kante,t)
    #-- Parameter
    (; R,P_max,U_soll,alpha,Jac_init) = kante.Param
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

    i_max = (UL-UR)/R; i_max = max(i_max,0);
    if haskey(Z,"P_max")==true
        i_max = min(i_max,P_max/UL);
    else
        i_max = min(i_max,kante.Param.i_max);
    end
    diff = U_soll - UR;
    ib = min(max(i,0),i_max);
    #fp = 0.5*(sign(diff)+1); fm = 1-fp; bound =  fp.*(i_max-ib)+fm.*(ib-0);
    bound = ifxaorb(diff,i_max-ib,ib);
    res = diff/alpha*bound;
    fp = 0.5*(sign(i-i_max)+1); fm = 0.5*(1+sign(0-i));
    res = res + 0.1*(fp*(i_max-i)+fm*(0-i)); #-- ziehe i ins intervall [0,i_max]
    #res = res +0.1*(ifxaorb(y_e-i_max,i_max-y_e,0.0)+ifxaorb(-y_e,-y_e,0.0));
    dy[k] = io*res-(1-io)*i;
end