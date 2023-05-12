Base.@kwdef mutable struct iBZ_Param
    F=96485; R=8.3145; z=2; DhO = 241.83e3; k=1.38e-23; 
    h=6.626e-34; Pstd = 101325; c=0.08988; # kg/m^3
    Eoc=65.7; Vl=58.4; Inom=8.128; Vnom=50.28; Imax=14.155; Vmin=45.707;
    N=65; #-- Anzahl Zellen
    nnom=58.83; #-- Effizienz
    Tc=42.3; Tnom=Tc+273; Pfuel=1.35; Pair=1.25; Vluftnom=14.91;
    xnom=0.9995; ynom=0.21; 
    NAnom=((Vl-Vnom)*(Imax-1)-(Vl-Vmin)*(Inom-1))/((log(Inom)*(Imax-1))-(log(Imax)*(Inom-1)));
    Rohm=(Vl-Vnom-NAnom*log(Inom))/(Inom-1); iOnom=exp((Vl-Eoc+Rohm)/(NAnom));
    Anom = NAnom /N; alpha = (R*Tnom)/(z*Anom*F);
    Uf_H2 =nnom*DhO*N/(100*(z*F*Vnom))  ;  # Gleichung 2-29
    Uf_02 =(60000*R*Tnom*N*Inom)/(4*F*Pair*Pstd*Vluftnom*ynom) ;# Gleichung 2-16 Pstd(normaler Luftdruck dazu für absoluten Druck
    PH2 = xnom*(1-Uf_H2)*Pfuel; P02 = ynom*(1-Uf_02)*Pair;
    Enomin = 1.229+(Tnom-298.15)*(-44.43/(z*F))+(R*Tnom/(z*F))*log(PH2*sqrt(P02));# Nominal für andere Werte
    Kl=2*F*k*(PH2*Pstd+P02*Pstd)/(h*R); Dg= -R*Tnom*log(iOnom/Kl);
    Ki=Eoc/Enomin; NAnomN=R*Tnom/(z*alpha*F)*N;  
    C1 = N*R*c; C2 = z*F*Uf_H2*xnom*Pstd;  C3 = z*F*k*Pstd/(R*h)*exp(-Dg/(R*Tnom)); 
    C4 = 1.229+(Tnom-298.15)*(-44.43)/(z*F); C5 = R*Tnom/(z*F);
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

    #-- Zustandsvariablen
    y = y_iBZ()

    #-- Stromknoten links und rechts
    KUL::Strom_Knoten
    KUR::Strom_Knoten

    #-- Gasknoten links und rechts
    KGL::Gas_Knoten
    KGR::Gas_Knoten

    #-- M-Matrix
    M::Array{Int} = [0; 0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::iBZ_kante,t)
    #-- Parameter
    (; C3,C4,C5,P02,Ki,NAnomN,Rohm,N,F,cv) = kante.Param
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
    P = P*1.0e-5; T = T_1 + 273.15; P = max(P,1.0e-5); 
    iOnomN = C3*(P+P02); EninReal = C4 + C5*log(max(P*sqrt(P02),1)); EocN = Ki*EninReal;

    #-- Kennlinie Brennstoffzelle
    #i0 = 1; s = 0.5*(1-cos(pi*min(i,i0)/i0)); s = min(s,1); #-- erst mal ohne Glättung 
    s=1;
    U_BZ = EocN - s*NAnomN*log(max(iBZ,iOnomN)/iOnomN)-Rohm*iBZ; #-- Spannung
    #U = max(U,0.0);
    #U = EocN - s*NAnomN*log(max(i/iOnomN,1.0e-1))-Rohm*i; #-- Spannung
    #m = C1*T*i/(C2*P); #-- Massenstrom, ???? Abhängigkeit vom Druck ???
    m_BZ = N*iBZ/(F*2)*2.01599e-3; #-- nach David

    dy[k] = io*(U_BZ-(UR-UL))+(1-io)*iBZ
    dy[k+1] = m - m_BZ
    dy[k+2] = e - cv*m*ifxaorb(m,TL,TR)
end
