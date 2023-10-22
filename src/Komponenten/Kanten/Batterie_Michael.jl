Base.@kwdef mutable struct iB_Param
    n_reih = 1
    n_par = 1
    A = 0.4919*n_reih;
    B = 130.2/3600/n_par;
   # koeff.K = 0.0224;
    K1 = 0.0224*n_reih;
    K2 = 0.0224/3600*n_reih;
    Q_max = 18.3*3600*n_par;
    R = 0.05*n_reih/n_par;
    U0 = 12.6481*n_reih;
   # alpha = koeff.K*koeff.Q_max/(3600*koeff.U0) 
    alpha = K2*Q_max/U0;
    soc_min = alpha/(1+alpha);
    SOC = 0.5
end

Base.@kwdef mutable struct y_iB
    Param::iB_Param
    i::Number = 0.0   
    u::Number = 0.0
    q::Number = Param.SOC*Param.Q_max
end

Base.@kwdef mutable struct iB_kante <: Strom_Kante
    #-- geänderte Parameter
    Param::iB_Param

    #-- Stromknoten links und rechts
    KL::Strom_Knoten
    KR::Strom_Knoten

    #-- Zustandsvariablen
    y = y_iB(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [0; 1; 1] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::iB_kante,t)
    #-- Parameter
    (; A,B,K1,K2,Q_max,R,U0,soc_min) = kante.Param

    #-- Zustandsvariablen
    iB = kante.y.i
    uB = kante.y.u
    qB = kante.y.q
    #--

    #-- Stromknoten links und rechts
    (; KL,KR) = kante
    UL = KL.y.U
    UR = KR.y.U
    #--
    
    uB = min(A,max(uB,0))
    soc = qB/Q_max; soc = minimum1(maximum1(soc,soc_min),1.095);
    U_B = U0 + uB - K1*ifxaorb(iB,1.0/soc,1.0/(1.1-soc))*iB - K2*(Q_max-qB)/soc;
    dy[k] = U_B - R*iB - (UR-UL);  #-- Bat.spannung 
    dy[k+1] = iB*B*(ifxaorb(iB,-uB,uB-A))
    dy[k+2] = -iB
    return nothing
end