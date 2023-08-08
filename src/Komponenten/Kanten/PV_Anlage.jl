Base.@kwdef mutable struct iPV_Param
    G_ref = 1000; 
    I_sc_stc = 3.11; I_ph_stc = I_sc_stc; 
    V_oc_stc = 21.8; Vmpp_stc = 17.0 ; Impp_stc = 2.88;
    a_sc = 0.0013; Tc_stc = 298.0; 
    eG = 1.12; # band energy gap in eV für Si
    Ns = 36.0; # Anzahl
    A_pv = 1.15;
    k_pv = 1.3805*1.0e-23; # Boltzmann - Konstante
    q = 1.6021*1.0e-19;
    faktor_a = A_pv*k_pv/q;
    Rs_pv = 0.45;  # müssen nochmal noch berechnet werden
    Rp = 310.0248;
    G = 2000.0
    T_PV = 298.15
    Module = 1
end

Base.@kwdef mutable struct y_iPV
    i::Number = 0.0
end

Base.@kwdef mutable struct iPV_kante <: Strom_Kante
    #-- geänderte Parameter
    Param::iPV_Param

    #-- Zustandsvariablen
    y = y_iPV()

    #-- Spannungsknoten links und rechts
    KL::Strom_Knoten
    KR::Strom_Knoten

    #-- M-Matrix
    M::Array{Int} = [0]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::iPV_kante,t)
    iPV = kante.y.i

    #-- Spannungsknoten links und rechts
    (; KL,KR,Z) = kante
    UL = KL.y.U
    UR = KR.y.U
    #--

    
    io = 1.0;  
    if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end
    # @show io # io Werte nochmal prüfen !!!  
    if (haskey(Z,"Leistung")==true)
        if isa(Z["Leistung"],Number) P = io * Z["Leistung"]; end #!!! io kommt nochmal in dy[k] = io*... vor!!!
        if isa(Z["Leistung"],Function) P = io * Z["Leistung"](t); end
        dy[k] = io*P/(UR-UL) - iPV
    elseif (haskey(Z,"Zeitreihe")==true)  
        wert = getwert(t,Z["zt"],Z["zwerte"],Z["interpol"],Z["ym"])
        t_scale = minimum1(t/60,1.0);
        P = wert*t_scale*kante.Param.Module
        dy[k] = io*P/(UR-UL) - iPV
    else #-- U-I Kennlinie
        (; G_ref,I_sc_stc,I_ph_stc,V_oc_stc,Vmpp_stc,Impp_stc,a_sc,Tc_stc,eG,Ns,faktor_a,Rs_pv,Rp,G,T_PV) = kante.Param

        U = UR - UL
        I_ph = G/G_ref*(I_ph_stc + a_sc*(T_PV - Tc_stc));
        a = Ns*faktor_a*T_PV;
        I0_ref = I_sc_stc/exp(V_oc_stc/a)
        I0 = I0_ref*(T_PV/Tc_stc)^3*exp(eG/faktor_a*(1/Tc_stc-1/T_PV))
        dy[k] = I_ph - I0*(exp((U+iPV*Rs_pv)/a)-1) - (U+Rs_pv*iPV)/Rp - iPV;
    end
end