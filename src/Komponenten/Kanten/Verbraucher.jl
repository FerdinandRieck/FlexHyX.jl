Base.@kwdef mutable struct iV_Param
    Scale = 1
end

Base.@kwdef mutable struct y_iV
    i::Number = 0.0
end

Base.@kwdef mutable struct iV_kante <: Strom_Kante
    Param::iV_Param

    #-- Spannungsknoten links und rechts
    KL::Strom_Knoten
    KR::Strom_Knoten

    #-- Zustandsvariablen
    y = y_iV()

    #-- M-Matrix
    M::Array{Int} = [0]

    #-- zusätzliche Infos
    Z::Dict
end

function Kante!(dy,k,kante::iV_kante,t)
    (; Scale) = kante.Param
    iV = kante.y.i
    Z = kante.Z

    #-- Spannungsknoten links und rechts
    (; KL,KR) = kante
    UL = KL.y.U
    UR = KR.y.U
    #--

    io = 1.0;  
    #!!! Wieso vorher => if get(Z,"Schaltzeit",0)!=0   ???
    if (haskey(Z,"Schaltzeit")==true) io = einaus(t,Z["Schaltzeit"],Z["Schaltdauer"]) end     
    # @show io # io Werte nochmal prüfen !!! 
    if (haskey(Z,"Leistung")==true)
        if isa(Z["Leistung"],Number) P = io * Z["Leistung"]; end #!!! io kommt nochmal in dy[k] = io*... vor!!!
        if isa(Z["Leistung"],Function) P = io * Z["Leistung"](t); end
    end
    if (haskey(Z,"Zeitreihe")==true)  
        wert = getwert(t,Z["zt"],Z["zwerte"],Z["interpol"],Z["ym"])
        P = wert*Scale
    end

    t_scale = min(1,t/1); #??? Wieso hier nicht t_scale = minimum1(t/60,1.0)??? Und wofür ist t_scale hier und normalerweise
    P = P*t_scale; 
    if (haskey(Z,"R")==true)
        dy[k] = iV - io*sqrt(abs(P)/Z["R"])*sign(P)
    else
        dy[k] = io*P/(UL-UR) - iV;
    end
end

function fcn_leistung(t)
    A = max(20.0,20.0+10*(t-8*3600)/(8*3600));
    P = 0.0+A*sin(t*2*pi/(4*3600))
    P = max(P,0.0)
end

function fcn_leistung_einaus(t)
	ts = range(3600.0,24*3600.0,step = 3600.0) 
	return 5.0*einaus(t,ts,10.0)
end
