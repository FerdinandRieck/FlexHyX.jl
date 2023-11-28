Base.@kwdef mutable struct E_Param
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
    Jac_init = true
end

Base.@kwdef mutable struct y_iE
    i::Number = 0.0
end

Base.@kwdef mutable struct y_mE
    m::Number = 0.0
    e::Number = 0.0
end

Base.@kwdef mutable struct iE_kante <: Strom_Kante
    #-- default Parameter
    Param::E_Param

    #-- Stromknoten links und rechts
    KL::Strom_Knoten
    KR::Strom_Knoten
 
    #-- mE_kante wird erst beim erzeugen der mE_kante übergeben
    K_mE = []

    #-- Zustandsvariablen
    y = y_iE()

    #-- M-Matrix
    M::Array{Int} = [0] 

    #-- zusätzliche Infos
    Z::Dict
end

Base.@kwdef mutable struct mE_kante <: Gas_Kante
    #-- default Parameter
    Param::E_Param

    #-- Gasknoten links und rechts
    KL::Gas_Knoten
    KR::Gas_Knoten 
    
    #-- iE_kante wird erst beim erzeugen übergeben
    K_iE = []

    #-- Zustandsvariablen
    y = y_mE()

    #-- M-Matrix
    M::Array{Int} = [0; 0] 

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::iE_kante,t)
    #-- Parameter
    (; A,r1,r2,s,t1,t2,t3,V_ref,n_Z,Jac_init) = kante.Param
    #--

    #-- Zustandsvariablen
    iE = kante.y.i;
    #--


    (; KL,KR,Z) = kante
    UL = KL.y.U
    UR = KR.y.U

    io = 1.0
    if Jac_init == true
        io = 1.0
    elseif (haskey(Z,"Schaltzeit")==true) 
        io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) 
    end

    T_L = 25 
    iA = iE/A 
    s1 = 0.5*(1-cos(pi*min(iE,1)))*s #-- Glättung
    U_el = (V_ref + (r1+r2*T_L)*iA + s1*log10((t1+t2/T_L+t3/(T_L^2))*max(iA,0)+1))*n_Z

    dy[k] = io*(U_el-(UL-UR))+(1-io)*iE
end

function Kante!(dy,k,kante::mE_kante,t)
    #-- Parameter
    (; nc,A,v_std,c,F,z,f_1m,f_1b,f_2b,f_2m,n_Z,cv,Jac_init) = kante.Param
    #--

    #-- Zustandsvariablen
    m = kante.y.m
    e = kante.y.e
    #--

    (; KL,KR,K_iE,Z) = kante
    TL = KL.y.T
    TR = KR.y.T
    iE = K_iE.y.i

    io = 1.0
    if Jac_init == true
        io = 1.0
    elseif (haskey(Z,"Schaltzeit")==true) 
        io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) 
    end

    T_L = 25 
    iA = iE/A 
    f_1 = f_1m*T_L + f_1b; f_2 = f_2b + f_2m*T_L;
    m_el = (v_std*c*iA^2/(f_1+iA^2)*f_2*iE*nc/(z*F))*n_Z;

    dy[k] = m - m_el
    dy[k+1] = e - cv*m*ifxaorb(m,TL,TR)
end

function iE_init!(knoten,kanten,M,kk,von,nach)
    kk["Typ"] = "E"
    Params = MakeParam(kk)
    kk["Typ"] = "iE"
    kante = iE_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk)
    push!(kanten, kante)
    append!(M, kante.M)
    if haskey(kk,"RefKante")
        error("“RefKante“ im JSON in “mE“ Kante eintragen!")
    end
end

function mE_init!(knoten,kanten,M,kk,von,nach)
    kk["Typ"] = "E"
    Params = MakeParam(kk)
    kk["Typ"] = "mE"
    kante = mE_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk)
    push!(kanten, kante)
    append!(M, kante.M)
    K_iE = kk["RefKante"]
    kanten[end].K_iE = kanten[K_iE]
    kanten[K_iE].K_mE = kanten[end]
end