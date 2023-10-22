Base.@kwdef mutable struct iBZ_Param
    T = 25.0 + 273.15;
    F=96485; R=8.3145; z=2; DhO = 241.83e3; k=1.38e-23; 
    h=6.626e-34; Pstd = 101325; c=0.08988; # kg/m^3
    # BZ
    N=36; #-- Anzahl Zellen
    nnom=46; #-- Effizienz
    Inom=51.8; Vnom=23.1;Tnom=25+273.15;Pfuel_nom=1.1; Pair_nom=1;Vluftnom=2400;
    xnom=0.9999; ynom=0.18;
    # Curvefit
    Anom = 0.04235;Eoc_nom =33;R_ohm=0.07568;i0_nom=1.008;
    
    alpha = (R*Tnom)/(z*Anom*F);
    Uf_H2 =nnom*DhO*N/(100*(z*F*Vnom));    # Gleichung 2-29
    Uf_02 =(60000*R*Tnom*N*Inom)/(4*F*Pair_nom*Pstd*Vluftnom*ynom); # Gleichung 2-16 Pstd(normaler Luftdruck dazu f?r absoluten Druck
    PH2_nom = xnom*(1-Uf_H2)*Pfuel_nom; P02_nom = ynom*(1-Uf_02)*Pair_nom;
    P02 = ynom*(1-Uf_02)*Pair_nom; # Luft ä®¤ert sich nicht
    E_nom = 1.229+(Tnom-298.15)*(-44.43/(z*F))+(R*Tnom/(z*F))*log(PH2_nom*sqrt(P02_nom));# Nominal f?r andere Werte
    Kl=2*F*k*(PH2_nom*Pstd+P02_nom*Pstd)/(h*R); Dg= -R*Tnom*log(i0_nom/Kl);
    Kc=Eoc_nom/E_nom; C1 = N*R*c; C2 = z*F*Uf_H2*xnom*Pstd;  C3 = z*F*k*Pstd/(R*h)*exp(-Dg/(R*T)); 
    C4 = 1.229+(T-298.15)*(-44.43)/(z*F); C5 = R*T/(z*F); C6=R/(z*alpha*F)*N;
    cv = 1.01798*1.0e4 
end

Base.@kwdef mutable struct y_iBZ
    i::Number = 0.0
    m::Number = 0.0
    e::Number = 0.0
end

Base.@kwdef mutable struct iBZ_kante <: Gas_Strom_Kante
    #-- default Parameter
    Param::iBZ_Param

    #-- Stromknoten links und rechts
    KUL::Strom_Knoten
    KUR::Strom_Knoten

    #-- Gasknoten links und rechts
    KGL::Gas_Knoten
    KGR::Gas_Knoten

    #-- Zustandsvariablen
    y = y_iBZ()

    #-- M-Matrix
    M::Array{Int} = [0; 0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::iBZ_kante,t)
    #-- Parameter
    (; xnom,Uf_H2,C3,C4,C5,C6,P02,Kc,R_ohm,N,F,cv) = kante.Param
    #--

    #-- Zustandsvariablen
    iBZ = kante.y.i;
    m = kante.y.m
    e = kante.y.e
    #--


    (; KUL,KUR,KGL,KGR,Z) = kante
    UL = KUL.y.U
    UR = KUR.y.U
    TL = KGL.y.T
    TR = KGR.y.T
    PL = KGL.y.P
    PR = KGR.y.P

    io = 1.0; if get(Z,"Schaltzeit",0)!=0 io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end

    P = PL-PR; T_1 = 25
    P = P*1.0e-5; P = max(P,1.0e-5); 
    T = T_1 + 273.15;
    PH2 = xnom*(1-Uf_H2)*P;  # schmiert bei 1.04 bar ab
    i0 = C3*(PH2+P02); En = C4 + C5*log(PH2*sqrt(P02)); Eoc = Kc*En;

    #-- Kennlinie Brennstoffzelle
    #-- i0 = 1; s = 0.5*(1-cos(pi*min(i,i0)/i0)); s = min(s,1); #-- erst mal ohne Glä´´ung 
    s=1
    NA=C6*T;
    U = Eoc - s*NA*log(max(iBZ*io,i0)/i0)-R_ohm*iBZ*io; #-- Spannung
    #m = C1*T*i/(C2*P); #-- Massenstrom, ???? Abhä®§igkeit vom Druck ???
    m_BZ = N*iBZ*io/(F*2)*2.01599e-3; #-- nach David

    #y_nR(i_ele) = U;
    KUR.y.U = U


    dy[k] = U-(UR-UL)
    dy[k+1] = m - m_BZ
    dy[k+2] = e - cv*m*ifxaorb(m,TL,TR)
end