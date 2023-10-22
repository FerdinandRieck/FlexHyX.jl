Base.@kwdef mutable struct iE_Param
    nc = 2;     #-- Anzahl Zellen in Parallelschaltung
    A=0.01; r1=3.538550e-4; r2=-3.02150e-6; s=0.22396;
    t1= 5.13093; t2=-2.40447e2; t3=3.410251e3;
    V_ref= 1.229;   #-- umkehrbare Spannung bei Standardbedinugnen
    v_std=0.0224136;    #-- IdealesGasvolumen bei Standardbedingungen
    c=0.08988;  #-- kg/m^3
    F=96485.33;     #-- Farraday- Konstante
    z=2;    #-- anzugebende Elektronen Wasser
    f_1m = 2.5; f_1b = 50; f_2b = 1; f_2m = -6.25e-6;
    n_Z = 1   #-- Anzahl Zellen  
    cv = 1.01798*1.0e4 
end

Base.@kwdef mutable struct y_iE
    i::Number = 0.0
    m::Number = 0.0
    e::Number = 0.0
end

Base.@kwdef mutable struct iE_kante <: Gas_Strom_Kante
    #-- default Parameter
    Param::iE_Param

    #-- Stromknoten links und rechts
    KUL::Strom_Knoten
    KUR::Strom_Knoten

    #-- Gasknoten links und rechts
    KGL::Gas_Knoten
    KGR::Gas_Knoten    

    #-- Zustandsvariablen
    y = y_iE()

    #-- M-Matrix
    M::Array{Int} = [0; 0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::iE_kante,t)
    #-- Parameter
    (; nc,A,r1,r2,s,t1,t2,t3,V_ref,v_std,c,F,z,f_1m,f_1b,f_2b,f_2m,n_Z,cv) = kante.Param
    #--

    #-- Zustandsvariablen
    iE = kante.y.i;
    m = kante.y.m
    e = kante.y.e
    #--


    (; KUL,KUR,KGL,KGR,Z) = kante
    UL = KUL.y.U
    UR = KUR.y.U
    TL = KGL.y.T
    TR = KGR.y.T

    io = 1.0; if get(Z,"Schaltzeit",0)!=0 io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end

    T_L = 25 #???Was ist T_L???
    iA = iE/A #-- ggf. ist iE = io*y_e besser?
    s1 = 0.5*(1-cos(pi*min(iE,1)))*s #-- Glättung
    U_el = (V_ref + (r1+r2*T_L)*iA + s1*log10((t1+t2/T_L+t3/(T_L^2))*max(iA,0)+1))*n_Z
    f_1 = f_1m*T_L + f_1b; f_2 = f_2b + f_2m*T_L;
    m_el = (v_std*c*iA^2/(f_1+iA^2)*f_2*iE*nc/(z*F))*n_Z;

    dy[k] = io*(U_el-(UL-UR))+(1-io)*iE
    dy[k+1] = m - m_el
    dy[k+2] = e - cv*m*ifxaorb(m,TL,TR)
end
