Base.@kwdef mutable struct GPMH_Param
    V_MH = 0.006; #-- Gesamtvolumen MHS in m^3
    epsilon = 0.44; #-- >Porosität
    V_s = (1-epsilon)*V_MH;  V_g = epsilon*V_MH;
    rho_max = 8300; rho_min = 8200; # kg/m^3
    delta_rho = rho_max- rho_min;
    M_leer = rho_min*V_s;
    delta_H = -30800; #Enthalpie [J/mol] 
    delta_S = -108; #Entropie [J/(mol K)]
    C_a = 2800; #reaction coefficient [1/s] [Mario]
    C_d = 9.57; 
    E_a = 31000; #activation energy [J/mol] [Mario] 
    E_d = 16450;
    MH2 = 0.00201588; #-- Molmasse von Wasserstoff in kg/mol
    R = 8.314; #-- universelle Gaskonstante in J/(mol K)
    Rs = R/MH2;
    p_A = 101325; #-- Atmosphärendruck in Pascal
    cp_s = 419; #-- J/(mol K) (Sophia)
    #cpH2 = 14890; %-- J/(mol K) (Sophia)
    cp_g = 14304; #-- spez. Wärmekapazitaet, nach David
    cv_g = cp_g - Rs; 
    L = 0.5; A = V_MH/L;  #-- Querschnitt, Länge
    K = 1.0e-8; mu = 8.813e-6; #-- Permeabilität, dyn. Viskosität
    KAepsdivmuL = K*A*epsilon/(mu*L/2); #-- für Massenfluss
    cv = 1.01798*1.0e4 
    P0 = 144210
    T0 = 293.15
    Θ0 = 0.5 
    rho_s = delta_rho * Θ0         
    rho_g = P0/(Rs*T0); 
    Masse0 = rho_g*V_g + rho_s*V_s;
end

Base.@kwdef mutable struct y_GPMH
    Param::GPMH_Param 
    M_H2::Number = Param.Masse0
    MT_MH::Number = (M_H2 + Param.M_leer)*Param.T0
    P::Number = Param.P0
    T::Number = Param.T0
    Θ::Number = Param.Θ0
end

Base.@kwdef mutable struct GPMH_Knoten <: Gas_Knoten
    #-- neue Parameter
    Param::GPMH_Param

    #-- Zustandsvariablen 
    y = y_GPMH(Param=Param)

    #-- M-Matrix
    M::Array{Int} = [1; 1; 0; 0; 1] 

    #-- zusätzeliche Infos
    Z::Dict
end

function Knoten!(dy,k,sum_i,sum_m,sum_e,knoten::GPMH_Knoten)
    #-- Parameter
    (; delta_rho,V_s,V_g,p_A,delta_H,R,delta_S,C_a,E_a,C_d,E_d,M_leer,cp_s,cv_g,MH2,Rs) = knoten.Param
    #--
    (; M_H2,MT_MH,P,T,Θ) = knoten.y
    
    rho_s = delta_rho*Θ; rho_g = (M_H2 - rho_s*V_s)/V_g
    p_eq = p_A*exp(delta_H/(R*T) - delta_S/R); #-- van't Hoff Gleichung
    A = C_a*exp(-E_a/(R*T))*log(max(P/p_eq,1.0))*(1-Θ);
    B = C_d*exp(-E_d/(R*T))*min(P-p_eq,0)/p_eq*Θ;
    dTheta_dt = ifxaorb(P-p_eq,A,B)
    cp = ((M_leer+V_s*rho_s)*cp_s + V_g*rho_g*cv_g)/(M_H2 + M_leer); #-- Gesamtwärmekapazität des MHS
    e_R = delta_rho*V_s*delta_H/MH2*dTheta_dt;

    dy[k] = sum_m
    dy[k+1] = (sum_e - e_R)/cp
    dy[k+2] = P-Rs*rho_g*T  
    dy[k+3] = T-MT_MH/(M_H2 + M_leer)
    dy[k+4] = dTheta_dt
end