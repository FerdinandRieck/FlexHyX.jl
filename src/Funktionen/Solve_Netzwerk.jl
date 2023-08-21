function solveNetzwerk(dir::String)
    println("---------------- This is FlexhyX ------------------")
#-- Netwerk einlesen
    J_cfg = JSON.parsefile(dir*"/Netzwerk/flexhyx.cfg")
    now = Dates.now(); jetzt = [Dates.year(now) Dates.month(now) Dates.day(now) Dates.hour(now) Dates.minute(now) 0]
    startzeit = get(J_cfg,"Startzeit",jetzt)
    startzeit = String(Symbol(startzeit'))
    startzeit = startzeit[5:end-1]
    startzeit = DateTime(startzeit,"yyyy mm dd HH MM SS")
    simdauer = get(J_cfg,"Simulationsdauer",86400)
    println("Startzeit, Simdauer: ",startzeit," ",simdauer)
    rtol = get(J_cfg,"RTOL",5.0e-4); atol = get(J_cfg,"ATOL",5.0e-4)
    println("rtol,atol: ",rtol," ",atol)
    pfad = get(J_cfg,"Pfad","."); pfad=dir*"/"*pfad*"/";
    netzfile = pfad*get(J_cfg,"Netzwerkfile",0)
    zeitfile = get(J_cfg,"Zeitreihenfile",nothing) 
    dtmax = get(J_cfg,"dtmax",600)

    znamen = []; zwerte = []; zt = [];
    if zeitfile != nothing
        zstart, zt, zwerte, znamen, zeinheit, ztitel = readZeitreihe(pfad*zeitfile)
        dt = Second(zstart-startzeit); dt = dt.value
        zt = zt .+ dt
    end

    knoten_infos, kanten_infos, eventfile = readNetz(dir, netzfile, zwerte, zt, znamen)

    n_n = size(knoten_infos)[1]; n_e = size(kanten_infos)[1];  
  
    M = Int[]; 
    kanten = Array{Any}(undef, n_e); knoten =  Array{Any}(undef, n_n);

    println("---------------- geänderte Parameter ------------------")

    for i = 1:n_n  #-- Knoten erzeugen ----------------------------
        kk = knoten_infos[i];  typ = kk["Typ"]; 

        #-- Parameter erzeugen und ändern
        Params = MakeParam(kk)
        #-- Knoten erzeugen
        s = Symbol(typ,"_Knoten"); obj = getfield(FlexHyX, s)
        knoten[i] = obj(Param=Params, Z=kk)     #-- z.B. U0_Knoten()

        M = append!(M, knoten[i].M)
    end

    #-- PW_max und TW_max suchen--------------------------
    PW_max = 0.0; TW_max = 0.0

    for i in eachindex(knoten_infos)
        if typeof(knoten[i]) <: Wasser_Knoten
            if hasfield(typeof(knoten[i].y), :P) == true
                PW_max = maximum([PW_max; knoten[i].y.P]) 
                TW_max = maximum([TW_max; knoten[i].y.T])
            end
        end
    end
    for i in eachindex(knoten_infos) #--- AW ändern ----
        kk = knoten[i].Z; typ = kk["Typ"];
        if typ=="WP" knoten[i].y.P = PW_max; end
        if typ=="T" knoten[i].y.T = TW_max; end
        # wieso kein T_max suchen? Wie können AW für Kopplungsknoten mit T besser bestimmt werden?
    end

    for i = 1:n_e  #-- Kanten erzeugen ---------------------------- 
        kk = kanten_infos[i]; typ = kk["Typ"]; 
        von = kk["VonNach"][1]; nach = kk["VonNach"][2]
        
        if typ=="iE"  #-- Elektrolyseeinheit
            mE = kk["RefKante"]; 
            von_mE = mE["VonNach"][1]; nach_mE = mE["VonNach"][2]
            Params = MakeParam(kk)
            kanten[i] = iE_kante(Param=Params, KUL=knoten[von], KUR=knoten[nach], KGL=knoten[von_mE], KGR=knoten[nach_mE], Z=kk)
        elseif typ=="iBZ"  #-- Brennstoffzelle
            mBZ = kk["RefKante"]; 
            von_mBZ = mBZ["VonNach"][1]; nach_mBZ = mBZ["VonNach"][2]
            Params = MakeParam(kk)
            kanten[i] = iBZ_kante(Param=Params, KUL=knoten[von], KUR=knoten[nach], KGL=knoten[von_mBZ], KGR=knoten[nach_mBZ], Z=kk)
        elseif typ=="mWTaR"  #-- Wärmetauscher
            Params = MakeParam(kk) 
            kanten[i] = mWTaR_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk) 
            if haskey(kk,"RohrAussen")
                RA = kk["RohrAussen"]
                kanten[i].RA = kanten[RA]
                kanten[RA].RA = kanten[i]
            end
        elseif typ=="mWTaR2"  #-- Wärmetauscher_reduziert
            Params = MakeParam(kk) 
            kanten[i] = mWTaR2_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk) 
            if haskey(kk,"RohrAussen")
                RA = kk["RohrAussen"]
                kanten[i].RA = kanten[RA]
                kanten[RA].RA = kanten[i]
            end
        elseif typ=="mWTaRK"  #-- Wärmetauscher_mdot_konstant
            Params = MakeParam(kk) 
            kanten[i] = mWTaRK_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk) 
            if haskey(kk,"RohrAustausch")
                R = kk["RohrAustausch"]
                kanten[i].R = kanten[R]
                kanten[R].R = kanten[i]
            end
        elseif typ=="mWTRK"  #-- Wärmeübertragung_mdot_konstant
            Params = MakeParam(kk) 
            kanten[i] = mWTRK_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk) 
            if haskey(kk,"KnotenAussen")
                KA = kk["KnotenAussen"]
                kanten[i].KA = knoten[KA]
                knoten[KA].in = kanten[[i]]
            end
        elseif typ=="mWTR"  #-- Wärmeübertragung_Rohr_reduziert
            Params = MakeParam(kk) 
            kanten[i] = mWTR_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk) 
            if haskey(kk,"KnotenAussen")
                KA = kk["KnotenAussen"]
                kanten[i].KA = knoten[KA]
                knoten[KA].in = kanten[[i]]
            end
        elseif typ=="mPI"  #-- PI-Regler
            Params = MakeParam(kk) 
            kanten[i] = mPI_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk) 
            K = kk["KnotenAussen"]
            kanten[i].K = knoten[K]
        elseif typ=="mPID"  #-- PID-Regler
            Params = MakeParam(kk) 
            kanten[i] = mPID_kante(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk) 
            K = kk["KnotenAussen"]
            kanten[i].K = knoten[K]
        else
            #-- Parameter erzeugen und ändern
            Params = MakeParam(kk) 
            #-- Kante erzeugen
            s = Symbol(typ,"_kante"); obj = getfield(FlexHyX, s)
            kanten[i] = obj(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk)    #-- z.B. iB_kante()
        end
        M = append!(M, kanten[i].M)
    end

    println("-------------------------------------------------------")

    #-- U_max, P_max suchen--------------------------
    U_max = 1.0; P_max = 101325;

    for i in eachindex(knoten)
        if hasfield(typeof(knoten[i].y), :U) == true
            U_max = maximum([U_max; knoten[i].y.U])
        end
        if (typeof(knoten[i]) <: Wasser_Knoten) == false
            if hasfield(typeof(knoten[i].y), :P) == true
                P_max = maximum([P_max; knoten[i].y.P]) 
            end
        end
    end
    for i in eachindex(kanten)
        if hasfield(typeof(kanten[i].y), :U) == true
            U_max = maximum([U_max; kanten[i].y.U]) 
        end
        if hasfield(typeof(kanten[i].Param), :U0) == true   #-- z. B. wegen U0 der Batterie
            U_max = maximum([U_max; kanten[i].Param.U0]) 
        end
        if (typeof(kanten[i]) <: Wasser_Kante) == false
            if hasfield(typeof(kanten[i].y), :P) == true
                P_max = maximum([P_max; kanten[i].y.P]) 
            end
        end
    end
    for i in eachindex(knoten) #--- AW ändern ----
        kk = knoten[i].Z; typ = kk["Typ"];
        if typ=="U" knoten[i].y.U = U_max; end 
        if typ=="GP" knoten[i].y.P = P_max; end 
        # wieso kein T_max suchen? Wie können AW für Kopplungsknoten mit T besser bestimmt werden?
    end
    #-----------------------------------------------------

    M = sparse(diagm(M))

    #-- Erzeuge Inzidenzmatrix 
    IM, IP = inzidenz(knoten,kanten) 
    
    for i in eachindex(knoten)
        in = findall(!iszero, IP[i,:])
        out = findall(!iszero, IM[i,:])
        if isempty(in) == false
            in
            knoten[i].in = kanten[in]
        end
        if isempty(out) == false
            knoten[i].out = kanten[out]    
        end
    end

    #-- Erzeuge Zustandsvektor y und Indizes wo was steht in y
    y, idx_iflussL, idx_iflussR, idx_mflussL, idx_mflussR, idx_eflussL, idx_eflussR, P_scale, idx_ele = netzwerk2array(knoten,kanten) 

    params = IM, IP, knoten, kanten, idx_iflussL, idx_iflussR, idx_mflussL, idx_mflussR, idx_eflussL, idx_eflussR, idx_ele

    println("Anzahl Gleichungen: ",length(y))

    #-- konsistente AW berechnen -----------
    ind_alg = findall(x->x==0,M[diagind(M)]);
    dy = 0*y;
    dgl!(dy,y,params,0.0);
    println("Test Vorher: ",Base.maximum(abs.(dy[ind_alg])))
    y_alg = copy(y[ind_alg])
    g!(dy_alg,y_alg) = f_aw!(dy_alg,y_alg,ind_alg,y,params)
    res = nlsolve(g!,y_alg)
    y[ind_alg] = res.zero;
    dgl!(dy,y,params,0.0);
    println("Test Nachher: ",Base.maximum(abs.(dy[ind_alg])))

    #--------------
    #-- Jacobi Struktur
    #f!(x,z) = dgl!(x,z,params,0.0); 
    #jac_sparsity = Symbolics.jacobian_sparsity(f!, similar(y), similar(y))  #-- funktioniert nicht immer
    dy0 = copy(y)
    jac_sparsity = Symbolics.jacobian_sparsity((dy, y) -> dgl!(dy, y, params, 0.0),dy0, y)
    #--------------

    t0 = time()
    f = ODEFunction(dgl!; mass_matrix=M)#, jac_prototype = float.(jac_sparsity))
    tspan = (0.0,simdauer)
    prob_ode = ODEProblem(f,y,tspan,params)

    if isempty(eventfile) 
        n_events = 0 
        sol = solve(prob_ode,Rodas5P(autodiff=true,diff_type=Val{:forward}),progress=true, reltol=rtol,abstol=atol,dtmax=dtmax)
    else
        global n_events
        cb = VectorContinuousCallback(event_condition,event_affect!,n_events,affect_neg! = nothing)
        sol = solve(prob_ode,Rodas5P(autodiff=true,diff_type=Val{:forward}), callback=cb, dense=false, progress=true, reltol=rtol, abstol=atol, dtmax=dtmax)
    end
    t1 = time()-t0
    println("CPU: ",t1)
    println(sol.retcode," nt=",size(sol.t)); 
    println(sol.destats)

    y = Leitsung_anhängen(sol,knoten,kanten,idx_iflussL,idx_iflussR,IM,IP)


    println("---------------- This was FlexHyX -----------------")
    return (idx_ele, sol, y, knoten_infos, kanten_infos, jac_sparsity)
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
    par = Dict("param" => D)
    Param = P(;(Symbol(k) => v for (k,v) in par["param"])...) #-- erstelle neue Param mit Änderungen
    return Param
end

function Leitsung_anhängen(y,knoten,kanten,idx_iflussL,idx_iflussR,IM,IP)

    idx2netzwerk!(knoten,kanten)
    
    y = Array(y)

    sum_i = IP[:,idx_iflussR[:,1]]*y[idx_iflussR[:,2],:] - IM[:,idx_iflussL[:,1]]*y[idx_iflussL[:,2],:];

    for i in eachindex(knoten)
        if (haskey(knoten[i].Z,"U0")==true)
            if (typeof(knoten[i]) == U0_Knoten) && (knoten[i].Z["U0"] > 0)
                idx_U = knoten[i].y.U
                U = y[[idx_U],:]
                I = sum_i[[i],:]
                P = U.*I
                y = vcat(y,P)
            end
        end
    end
    for i in eachindex(kanten)
        if (typeof(kanten[i]) <: Strom_Kante) && (typeof(kanten[i]) != iSP0_kante)
            idx_UL = kanten[i].KL.y.U
            idx_UR = kanten[i].KR.y.U
            idx_i  = kanten[i].y.i
            UL = y[[idx_UL],:]; UR = y[[idx_UR],:]
            U = UR - UL; I = y[[idx_i],:]
            P = U.*I
            y = vcat(y,P)
        elseif typeof(kanten[i]) == iSP0_kante
            idx_UL = kanten[i].KL.y.U
            idx_UR = kanten[i].KR.y.U
            idx_iL  = kanten[i].y.iL
            idx_iR  = kanten[i].y.iR
            UL = y[[idx_UL],:]; UR = y[[idx_UR],:]
            IL = y[[idx_iL],:]; IR = y[[idx_iR],:]
            P = UR.*IR - UL.*IL;
            y = vcat(y,P)
        elseif typeof(kanten[i]) <: Gas_Strom_Kante
            idx_UL = kanten[i].KUL.y.U
            idx_UR = kanten[i].KUR.y.U
            idx_i  = kanten[i].y.i
            UL = y[[idx_UL],:]; UR = y[[idx_UR],:]
            U = UR - UL; I = y[[idx_i],:]
            P = U.*I
            y = vcat(y,P)
        end
    end
    return y
end


