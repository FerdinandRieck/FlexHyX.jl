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
    println("Startzeit, Simdauer:",startzeit," ",simdauer)
    rtol = get(J_cfg,"RTOL",5.0e-4); atol = get(J_cfg,"ATOL",5.0e-4)
    println("rtol,atol:",rtol," ",atol)
    pfad = get(J_cfg,"Pfad","."); pfad=dir*"/"*pfad*"/";
    netzfile = pfad*get(J_cfg,"Netzwerkfile",0)
    zeitfile = get(J_cfg,"Zeitreihenfile",nothing) 

    znamen = []; zwerte = []; zt = [];
    if zeitfile != nothing
        zstart, zt, zwerte, znamen, zeinheit, ztitel = readZeitreihe(pfad*zeitfile)
        dt = Second(zstart-startzeit); dt = dt.value
        zt = zt .+ dt
    end

    knoten_infos, kanten_infos, eventfile = readNetz(dir, netzfile, zwerte, zt, znamen)

    #-- Anfangswerte setzen
    IM, IP = inzidenz(knoten_infos,kanten_infos)

    n_n = size(knoten_infos)[1]; n_e = size(kanten_infos)[1];  
  
    M = Int[]; 
    kanten = Array{Any}(undef, n_e); knoten =  Array{Any}(undef, n_n); 

    println("---------------- geänderte Parameter ------------------")

    for i=1:n_n  #-- Knoten erzeugen ----------------------------
        kk = knoten_infos[i];  typ = kk["Typ"]; 

        #-- Parameter erzeugen und ändern
        Params = MakeParam(kk)
        #-- Knoten erzeugen
        s = Symbol(typ,"_Knoten"); obj = getfield(FlexHyX, s)
        knoten[i] = obj(Param=Params, Z=kk)     #-- z.B. U0_Knoten()

        M = vcat(M, knoten[i].M)
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
        else
            #-- Parameter erzeugen und ändern
            Params = MakeParam(kk) 
            #-- Kante erzeugen
            s = Symbol(typ,"_kante"); obj = getfield(FlexHyX, s)
            kanten[i] = obj(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk)    #-- z.B. iB_kante()
        end

        M = vcat(M, kanten[i].M)
    end

    println("-------------------------------------------------------")

    #-- U_max, P_max, PW_max suchen--------------------------
    U_max = 0; P_max = 0; PW_max = 0

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
        if typ=="U" knoten[i].y.U = U_max; end 
        if typ=="GP" knoten[i].y.P = P_max; end 
        if typ=="WP" knoten[i].y.PW = PW_max; end
    end
    #-----------------------------------------------------

    M = sparse(diagm(M))

    #-- Erzeuge Zustandsvektor y und Indizes wo was steht in y 
    elemente = Netzwerk(kanten=kanten,knoten=knoten)  #-- gesamtes Netzwerk  

    y, idx_iflussL, idx_iflussR, idx_mfluss, idx_efluss, P_scale, leg_knoten, leg_kanten, idx_ele = netzwerk2array(elemente)  

    params = IM, IP, elemente, idx_iflussL, idx_iflussR, idx_mfluss, idx_efluss, idx_ele, n_n, n_e
             
    #-- konsistente AW berechnen -----------
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
        sol = solve(prob_ode,Rodas5P(autodiff=true,diff_type=Val{:forward}),progress=true, reltol=rtol,abstol=atol,dtmax=600)
    else
        global n_events
        cb = VectorContinuousCallback(event_condition,event_affect!,n_events,affect_neg! = nothing)
        sol = solve(prob_ode,Rodas5P(autodiff=true,diff_type=Val{:forward}), callback=cb, dense=false, progress=true, reltol=rtol, abstol=atol, dtmax=600)
    end

    y = Leitsung_anhängen(sol,elemente,idx_iflussL,idx_iflussR,IM,IP)

    t1 = time()-t0
    println("CPU:",t1)
    println(sol.retcode," nt=",size(sol.t)); 
    println(sol.destats)
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

function Leitsung_anhängen(y,elemente,idx_iflussL,idx_iflussR,IM,IP)

    idx2netzwerk!(elemente)
    
    y = Array(y)

    sum_i = IP[:,idx_iflussR[:,1]]*y[idx_iflussR[:,2],:] - IM[:,idx_iflussL[:,1]]*y[idx_iflussL[:,2],:];

    for i = 1:length(elemente.knoten)
        if (haskey(elemente.knoten[i].Z,"U0")==true)
            if (typeof(elemente.knoten[i]) == U0_Knoten) && (elemente.knoten[i].Z["U0"] > 0)
                idx_U = elemente.knoten[i].y.U
                U = y[[idx_U],:]
                I = sum_i[[i],:]
                P = U.*I
                y = vcat(y,P)
            end
        end
    end
    for i = 1:length(elemente.kanten)
        if (typeof(elemente.kanten[i]) <: Strom_Kante) && (typeof(elemente.kanten[i]) != iSP0_kante)
            idx_UL = elemente.kanten[i].KL.y.U
            idx_UR = elemente.kanten[i].KR.y.U
            idx_i  = elemente.kanten[i].y.i
            UL = y[[idx_UL],:]; UR = y[[idx_UR],:]
            U = UR - UL; I = y[[idx_i],:]
            P = U.*I
            y = vcat(y,P)
        elseif typeof(elemente.kanten[i]) == iSP0_kante
            idx_UL = elemente.kanten[i].KL.y.U
            idx_UR = elemente.kanten[i].KR.y.U
            idx_iL  = elemente.kanten[i].y.i
            idx_iR  = elemente.kanten[i].y.i_out
            UL = y[[idx_UL],:]; UR = y[[idx_UR],:]
            IL = y[[idx_iL],:]; IR = y[[idx_iR],:]
            P = UR.*IR - UL.*IL;
            y = vcat(y,P)
        elseif typeof(elemente.kanten[i]) <: Gas_Strom_Kante
            idx_UL = elemente.kanten[i].KUL.y.U
            idx_UR = elemente.kanten[i].KUR.y.U
            idx_i  = elemente.kanten[i].y.i
            UL = y[[idx_UL],:]; UR = y[[idx_UR],:]
            U = UR - UL; I = y[[idx_i],:]
            P = U.*I
            y = vcat(y,P)
        end
    end
    return y
end