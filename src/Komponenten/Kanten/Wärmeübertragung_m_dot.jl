Base.@kwdef mutable struct mWT_Param
    #---Wärmeübertragung ----
    #-- Berechnung des Wärmedurchgangskoeffizienten (nach VDI-Wärmeatlas, 11.
    #-- Auflage, Kapitel G2)
    #-- Stoffwerte
    nu = 1e-6; # kinematische Viskosität von Wasser (bei 20°C) in m^2/s
    kappa = 1.4e-7; # Temperaturleitfähigkeit von Wasser (bei 20°C) in m^2/s
    lambda_w = 0.604; # Wärmeleitfähigkeit von Wasser (bei 20°C) in W/(m K)
    lambda_f = 236; # Wärmeleitfähigkeit des Feststoffs (Wand) in W/(m K), Wert von Sophia
    rho_w =998.21; # Dichte von Wasser (bei 20°C) in kg/m^3
    cp_w = 4181; # spezifische Wärmekapazität von Wasser (bei 20°C) in J/(kg K)
    #-- Stellgrößen (bisher willkürlich gesetzt)
    di = 0.3;   # Innendurchmesser des Wärmeübertragers in m
    da = 0.35;  # Außendurchmesser des Wärmeübertragers in m
    L = 1; # Länge des Wärmeübertragers in m
    delta = 3.76e-3; # Wandstärke in m, Wert von Sophia
    #--Berechnungen
    dh = da - di;   # hydraulischer Durchmesser
    A_q = pi/2*(da^2 - di^2);   # Querschnittsfläche des Wärmeübertragers
    Aa = da*pi*L;
    Ai = di*pi*L;
    Am = (Aa - Ai)/log(Aa/Ai);
    Pr = nu/kappa;
    Nu1 = 3.66 + 1.2*(di/da)^(-0.8);
    fg = 1.615*(1 + 0.14*(di/da)^(-0.5));
    #-- ab hier abhängig vom Massenstrom Wasser
    V_w = 0.0001; # Volumenstrom des Kühlwassers in m^3/s
    Re = V_w*dh/(A_q*nu);   # Reynolds-Zahl
    Nu2 = fg*(Re*Pr*dh/L)^(1/3);
    Nu3 = (2/(1+22*Pr))^(1/6)*(Re*Pr*dh/L)^(1/2);
    Nu = (Nu1^3 + Nu2^3 + Nu3^3)^(1/3);
    alpha = Nu*lambda_w/dh;
    #kA = 1/(delta/(lambda_f*Am) + 1/(alpha*Aa));
    mcp_w = V_w*rho_w*cp_w;
    #---
    # Q_dot = mcp_w*(TW_zu - T_MHS)*(1 - exp(-kA/mcp_w));
    kA = 5;
    # Alpha = Defaultwert festlegen
    # A = Defaultwert festlegen
end

Base.@kwdef mutable struct y_mWT
    m::Number = 0.0
    eL::Number = 0.0
    eR::Number = 0.0
end

Base.@kwdef mutable struct mWT_kante <: Temp_Kante
    #-- default Parameter
    Param::mWT_Param

    #-- Zustandsvariablen
    y = y_mWT()

    #-- Gasknoten links und rechts
    KL::Wasser_Knoten  
    KR::Wasser_Knoten  

    #-- Rohr links 
    RL = 0 

    #-- M-Matrix
    M::Array{Int} = [0; 0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::mWT_kante,t)
    #-- Parameter
    (; kA) = kante.Param
    #--

    #-- Zustandsvariable
    eL = kante.y.eL
    eR = kante.y.eR
    m = kante.y.m
    #--

    (; KL,KR,RL,Z) = kante
    TL = KL.y.T
    TR = KR.y.T
    m_in = RL.y.mR
    e_in = RL.y.eR

    T_aussen = 263.15


    dy[k] = m - m_in
    dy[k+1] = eL - e_in
    if (haskey(Z,"alpha")==true) && (haskey(Z,"A")==true)
        dy[k+2] = eR - Z["alpha"]*Z["A"]*(TL-TR)
    else
        dy[k+2] = eR - (e_in - kA*(T_aussen-TL))  #-- hier noch dot m einfügen
    end
end