Base.@kwdef mutable struct iB_Param
    A = 0.4919; 
    B = 1.302; 
    K1 = 0.0224; 
    K2 = 6.222222222222222e-06; 
    Q_max = 65880; 
    R = 0.5; 
    U0 = 12.6481; 
    soc_min = 0.031392201880530
    SOC = 0.5
end

Base.@kwdef mutable struct y_iB
    Param::iB_Param
    iB::Number = 0.0   
    uB::Number = 0.0
    qB::Number = Param.SOC*Param.Q_max
end

Base.@kwdef mutable struct iB_kante <: Strom_Kante
    #-- geänderte Parameter
    Param::iB_Param

    #-- Zustandsvariablen
    y = y_iB(Param=Param)

    #-- Spannungsknoten links und rechts
    KL::Strom_Knoten
    KR::Strom_Knoten

    #-- M-Matrix
    M::Array{Int} = [0; 1; 1] 

    #-- Jacobi Struktur
    J::Array{Int} = [1 1 1; 1 1 1; 1 1 1]
    J_KL::Dict = Dict("eq1" => ["U"], "eq2" => ["U"], "eq3" => ["U"]) #["P" "T"; "T" "0"; "0"]
    J_KR::Dict = Dict("eq1" => ["U"], "eq2" => ["U"], "eq3" => ["U"])
    
    #-- zusätzliche Infos
    Z::Dict

    #-- Ladestatus für Event-Funktionen eventuell auch hier hin statt in Z
    #status_laden::Int = 0
    #status_entladen::Int = 0
end

function Kante!(dy,k,kante::iB_kante,t)
    #-- Parameter
    (; A,B,K1,K2,Q_max,R,U0,soc_min) = kante.Param

    #-- Zustandsvariablen
    (; iB,uB,qB) = kante.y
    #--

    #-- Spannungsknoten links und rechts
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