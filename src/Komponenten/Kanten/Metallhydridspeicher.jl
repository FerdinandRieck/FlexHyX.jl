Base.@kwdef mutable struct mMH_Param
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

Base.@kwdef mutable struct y_mMH
    m::Number = 0
    e::Number = 0
end

Base.@kwdef mutable struct mMH_kante <: Gas_Kante
    #-- default Parameter
    Param::mMH_Param

    #-- Gasknoten links und rechts
    KL::Gas_Knoten
    KR::GPMH_Knoten

    #-- Zustandsvariablen
    y = y_mMH()

    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mMH_kante,t)
        #-- Parameter
        (; KAepsdivmuL, delta_rho,V_s,V_g,cv) = kante.Param
        #--
    
        #-- Zustandsvariablen 
        m = kante.y.m
        e = kante.y.e
        #--
    
        #-- Knoten links und rechts
        (; KL,KR) = kante
        PL = KL.y.P
        PR = KR.y.P
        TL = KL.y.T
        TR = KR.y.T
        Θ = KR.y.Θ
        M_H2 = KR.y.M_H2
        #--

        rho_s = delta_rho*Θ
        rho_g = (M_H2 - rho_s*V_s)/V_g
        delta_P = PL - PR;
        md = KAepsdivmuL*rho_g*delta_P;  #--- Massenfluss

        dy[k] = m - md
        dy[k+1] = e - cv*m*ifxaorb(m,TL,TR)
end
