Base.@kwdef mutable struct BZ_Param
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
    Jac_init = true
end

Base.@kwdef mutable struct y_iBZ
    i::Number = 0.0
end

Base.@kwdef mutable struct y_mBZ
    m::Number = 0.0
    e::Number = 0.0
end

Base.@kwdef mutable struct iBZ_kante <: Strom_Kante
    #-- default Parameter
    Param::BZ_Param

    #-- Stromknoten links und rechts
    KL::Strom_Knoten
    KR::Strom_Knoten

    #-- mBZ_kante wird erst beim erzeugen der mBZ_kante übergeben
    K_mBZ = []

    #-- Zustandsvariablen
    y = y_iBZ()

    #-- M-Matrix
    M::Array{Int} = [0] 

    #-- zusätzliche Infos
    Z::Dict
end

Base.@kwdef mutable struct mBZ_kante <: Gas_Kante
    #-- default Parameter
    Param::BZ_Param

    #-- Gasknoten links und rechts
    KL::Gas_Knoten
    KR::Gas_Knoten

    #-- iBZ_kante wird erst beim erzeugen übergeben
    K_iBZ = []

    #-- Zustandsvariablen
    y = y_mBZ()

    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::iBZ_kante,t)
    #-- Parameter
    (; xnom,Uf_H2,C3,C4,C5,C6,P02,Kc,R_ohm,Jac_init) = kante.Param
    #--

    #-- Zustandsvariablen
    iBZ = kante.y.i;
    #--

    (; KL,KR,K_mBZ,Z) = kante
    UL = KL.y.U
    UR = KR.y.U
    PL = K_mBZ.KL.y.P
    PR = K_mBZ.KR.y.P

    io = 1.0
    if Jac_init == true
        io = 1.0
    elseif (haskey(Z,"Schaltzeit")==true) 
        io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) 
    end

    P = PL-PR; T_1 = 25
    P = P*1.0e-5; P = max(P,1.0e-5); 
    T = T_1 + 273.15;
    PH2 = xnom*(1-Uf_H2)*P; 
    i0 = C3*(PH2+P02); En = C4 + C5*log(PH2*sqrt(P02)); Eoc = Kc*En;

    #-- Kennlinie Brennstoffzelle
    s=1
    NA=C6*T;
    U = Eoc - s*NA*log(max(iBZ*io,i0)/i0)-R_ohm*iBZ*io; #-- Spannung

    #Q_dot = max(0,min(1.0e3,1.12*(1.25 + U)*iBZ))

    dy[k] = U-(UR-UL)
end

function Kante!(dy,k,kante::mBZ_kante,t)
    #-- Parameter
    (; N,F,cv,Jac_init) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    e = kante.y.e
    #--

    (; KL,KR,K_iBZ,Z) = kante
    TL = KL.y.T
    TR = KR.y.T
    iBZ = K_iBZ.y.i

    io = 1.0
    if Jac_init == true
        io = 1.0
    elseif (haskey(Z,"Schaltzeit")==true) 
        io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) 
    end

    m_BZ = N*iBZ*io/(F*2)*2.01599e-3; #-- nach David

    #Q_dot = max(0,min(1.0e3,1.12*(1.25 + U)*iBZ))

    dy[k] = m - m_BZ
    dy[k+1] = e - (cv*m*ifxaorb(m,TL,TR))
end

function iBZ_init(knoten,kanten,M,kk,von,nach)
    kk["Typ"] = "BZ"
    Params = MakeParam(kk)
    kk["Typ"] = "iBZ"
    kante = iBZ_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk)
    push!(kanten,kante)
    append!(M, kante.M)
end

function mBZ_init(knoten,kanten,M,kk,von,nach)
    kk["Typ"] = "BZ"
    Params = MakeParam(kk)
    kk["Typ"] = "mBZ"
    kante = mBZ_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk)
    push!(kanten,kante)
    append!(M, kante.M)
    K_iBZ = kk["RefKante"]
    kanten[end].K_iBZ = kanten[K_iBZ]
    kanten[K_iBZ].K_mBZ = kanten[end]
end