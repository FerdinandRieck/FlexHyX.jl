function solveNetzwerk(dir::String)
    println("---------------- This is FlexhyX ------------------")
#-- Netwerk einlesen
    #dir = dirname(@__DIR__)
    J_cfg = JSON.parsefile(dir*"/Netzwerk/flexhyx.cfg")
    now = Dates.now(); jetzt = [Dates.year(now) Dates.month(now) Dates.day(now) Dates.hour(now) Dates.minute(now) 0]
    startzeit = get(J_cfg,"Startzeit",jetzt)
    startzeit = String(Symbol(startzeit'))
    startzeit = startzeit[5:end-1]
    startzeit = DateTime(startzeit,"yyyy mm dd HH MM SS")
    simdauer = get(J_cfg,"Simulationsdauer",86400)
    println("Startzeit, Simdauer:",startzeit," ",simdauer)
    rtol = get(J_cfg,"RTOL",5.0e-4); atol = get(J_cfg,"ATOL",5.0e-4)
    println("rtol,atol:",rtol," ",atol)
    pfad = get(J_cfg,"Pfad","."); pfad=dir*"/"*pfad*"/";
    netzfile = pfad*get(J_cfg,"Netzwerkfile",0)
    zeitfile = get(J_cfg,"Zeitreihenfile",nothing) 

    znamen = []; zwerte = []; zt = [];
    if zeitfile != nothing
        zstart, zt, zwerte, znamen, zeinheit, ztitel = read_zeitreihe(pfad*zeitfile)
        dt = Second(zstart-startzeit); dt = dt.value
        zt = zt .+ dt
    end

    eventfile, knoten_infos, kanten_infos = read_netz(dir,netzfile, zwerte, zt, znamen)

   
    #-- Anfangswerte setzen
    IM, IP = inzidenz(knoten_infos,kanten_infos)
    n_n = size(knoten_infos)[1]; n_e = size(kanten_infos)[1];  
  
    M = Int[]; 
    kanten = Array{Any}(undef, n_e); knoten =  Array{Any}(undef, n_n); 

    U_max = 0; P_max = 0; PW_max = 0

    println("---------------- geänderte Parameter ------------------")

    for i=1:n_n  #-- Knoten erzeugen ----------------------------
        kk = knoten_infos[i];  typ = kk["Typ"]; 

        #-- Parameter erzeugen und ändern
        Params = MakeParam(kk)
        #-- Knoten erzeugen
        s = Symbol(typ,"_Knoten"); obj = getfield(FlexHyX, s)
        knoten[i] = obj(Param=Params, Z=kk)     #-- z.B. U0_Knoten()

        M = [M; knoten[i].M]
    end

    for i=1:n_e  #-- Kanten erzeugen ---------------------------- 
        kk = kanten_infos[i]; typ = kk["Typ"]; 
        von = kk["VonNach"][1]; nach = kk["VonNach"][2]
        
        if typ=="iE"
            mE = kk["RefKante"]; 
            von_mE = mE["VonNach"][1]; nach_mE = mE["VonNach"][2]
            Params = MakeParam(kk)
            kanten[i] = iE_kante(Param=Params, KUL=knoten[von], KUR=knoten[nach], KGL=knoten[von_mE], KGR=knoten[nach_mE], Z=kk)
        elseif typ=="iBZ"
            mBZ = kk["RefKante"]; 
            von_mBZ = mBZ["VonNach"][1]; nach_mBZ = mBZ["VonNach"][2]
            Params = MakeParam(kk)
            kanten[i] = iBZ_kante(Param=Params, KUL=knoten[von], KUR=knoten[nach], KGL=knoten[von_mBZ], KGR=knoten[nach_mBZ], Z=kk)
        elseif typ=="mMH" 
            kk_GPMH = knoten[nach].Z    #-- Infos von GPMH_Knoten
            Params = MakeParam(kk_GPMH)
            knoten[nach] = GPMH_Knoten(Param=Params, KL=knoten[von], Z=kk_GPMH) # GPMH_Knoten neu anlegen mit neuen Params und KL
            kk["Typ"] = "GPMH"  #-- Damit selber Parametersatz wie von GPMH_Param verwendet wird
            Params = MakeParam(kk)
            kanten[i] = mMH_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk)
        else
            #-- Parameter erzeugen und ändern
            Params = MakeParam(kk) 
            #-- Kante erzeugen
            s = Symbol(typ,"_kante"); obj = getfield(FlexHyX, s)
            kanten[i] = obj(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk)    #-- z.B. iB_kante()
        end

        M = [M; kanten[i].M] 
    end

    println("-------------------------------------------------------")

    #-- U_max, P_max, PW_max suchen--------------------------
    for i=1:n_n
        if hasfield(typeof(knoten[i].y), :U) == true
            U_max = max(U_max,knoten[i].y.U) 
        end
        if hasfield(typeof(knoten[i].y), :P) == true
            P_max = max(P_max,knoten[i].y.P) 
        end
        if hasfield(typeof(knoten[i].y), :PW) == true
            PW_max = max(PW_max,knoten[i].y.PW) 
        end
    end
    for i=1:n_e
        if hasfield(typeof(kanten[i].y), :U) == true
            U_max = max(U_max,kanten[i].y.U) 
        end
        if hasfield(typeof(kanten[i].Param), :U0) == true   #-- z. B. wegen U0 der Batterie
            U_max = max(U_max,kanten[i].Param.U0) 
        end
        if hasfield(typeof(kanten[i].y), :P) == true
            P_max = max(P_max,kanten[i].y.P) 
        end
        if hasfield(typeof(kanten[i].y), :PW) == true
            PW_max = max(PW_max,kanten[i].y.PW) 
        end
    end
    for i=1:n_n #--- AW ändern ----
        kk = knoten[i].Z; typ = kk["Typ"];
        if typ=="U" knoten[i].y.U = U_max; end  #???Wofür wird das genau benötigt AW??? was wenn zwei nicht gekoppelte Stromnetze???
        if typ=="GP" knoten[i].y.P = P_max; end #???Wieso kein T_max bestimmen???
        if typ=="WP" knoten[i].y.PW = PW_max; end
    end
    #-----------------------------------------------------

    M = sparse(diagm(M))

    #-- Erzeuge Zustandsvektor y und Indizes wo was steht in y 
    elemente = Netzwerk(kanten=kanten,knoten=knoten)  #-- gesamtes Netzwerk 

    y, idx_iflussL, idx_iflussR, idx_mfluss, idx_efluss, P_scale, y_leg, idx_ele = tuple2array(elemente)  

   # Jacstru = comp_jacstru(IP,IM,idx_ifluss,idx_mfluss,idx_efluss,elemente,length(y))

    #-- Netzinfo und Speicherplatz (übergebe als Parameter an solver/dgl-function)
    i_flussL = Array{Number}(undef, n_e);  #-- nur einmal Speicher reservieren
    i_flussR = Array{Number}(undef, n_e);
    m_fluss = Array{Number}(undef, n_e); 
    e_fluss = Array{Number}(undef, n_e); 

    #-------------

    params = IM, IP, elemente, i_flussL, i_flussR, m_fluss, e_fluss, idx_iflussL, idx_iflussR, idx_mfluss, idx_efluss, idx_ele, n_n, n_e
             
    #--------------
    ind_alg = findall(x->x==0,M[diagind(M)]);
    dy = 0*y;
    dgl!(dy,y,params,0.0);
    println("Test Vorher:",Base.maximum(abs.(dy[ind_alg])))
    y_alg = copy(y[ind_alg])
    g!(dy_alg,y_alg) = f_aw!(dy_alg,y_alg,ind_alg,y,params)
    res = nlsolve(g!,y_alg)
    y[ind_alg] = res.zero;
    dgl!(dy,y,params,0.0);
    println("Test Nachher:",Base.maximum(abs.(dy[ind_alg])))

    #--------------
    t0 = time()
    f = ODEFunction(dgl!,mass_matrix=M)
    tspan = (0.0,simdauer)
    prob_ode = ODEProblem(f,y,tspan,params)

    if isempty(eventfile) 
        n_events = 0 
        sol = solve(prob_ode,Rodas5P(autodiff=true,diff_type=Val{:forward}),progress=true, reltol=rtol,abstol=atol,dtmax=600)
    else
        global n_events
        cb = VectorContinuousCallback(event_condition,event_affect!,n_events,affect_neg! = nothing)
        sol = solve(prob_ode,Rodas5P(autodiff=true,diff_type=Val{:forward}), callback=cb, dense=false, progress=true, reltol=rtol, abstol=atol, dtmax=600)
    end

    t1 = time()-t0
    println("CPU:",t1)
    println(sol.retcode," nt=",size(sol.t)); 
    println(sol.destats)
    println("---------------- This was FlexHyX -----------------")
    return sol
end

function MakeParam(kk) 
    s = Symbol(kk["Typ"],"_Param"); P = getfield(FlexHyX, s)
    Param = P() #-- erzeuge Param z.B. iB_Param
    D = Dict()
    for ff in fieldnames(typeof(Param))
        if haskey(kk,String(ff))==true 
            D[String(ff)] = kk[String(ff)] #-- speichere geänderte Params in Dict
            println(kk["Typ"]," --> ",String(ff),"=",kk[String(ff)])
        end
    end
    par = Dict("param"=>D)
    Param = P(;(Symbol(k) => v for (k,v) in par["param"])...) #-- erstelle neue Param mit Änderungen
    return Param
end