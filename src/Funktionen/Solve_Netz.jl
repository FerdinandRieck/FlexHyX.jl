function solveNetz(dir::String)
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
    AW = get(J_cfg,"AW",1)
    JPattern = get(J_cfg,"JPattern",1)

    znamen = []; zwerte = []; zt = [];
    if zeitfile != nothing
        zstart, zt, zwerte, znamen, zeinheit, ztitel = readZeitreihe(pfad*zeitfile)
        dt = Second(zstart-startzeit); dt = dt.value
        zt = zt .+ dt
    end

    knoten_infos, kanten_infos, eventfile = readNetz(dir, netzfile, zwerte, zt, znamen)
    n_n = size(knoten_infos)[1]; n_e = size(kanten_infos)[1];  
    M = Int[]; kanten = []; knoten = [];

    #println("---------------- geänderte Parameter ------------------")

    for i = 1:n_n  #-- Knoten erzeugen ----------------------------
        kk = knoten_infos[i];  typ = kk["Typ"]; 

        #-- Parameter erzeugen und ändern
        Params = MakeParam(kk)
        #-- Knoten erzeugen
        s = Symbol(typ,"_Knoten"); obj = getfield(FlexHyX, s)
        push!(knoten,obj(Param=Params, Z=kk))     #-- z.B. U0_Knoten()
        append!(M, knoten[end].M)
    end

    for i = 1:n_e  #-- Kanten erzeugen ---------------------------- 
        kk = kanten_infos[i]; typ = kk["Typ"]; 
        von = kk["VonNach"][1]; nach = kk["VonNach"][2]

        if isdefined(FlexHyX, Symbol(typ,"_init!"))
            #-- spezielle Kanten erzeugen
            s = Symbol(typ,"_init!"); obj_init! = getfield(FlexHyX, s)
            obj_init!(knoten,kanten,M,kk,von,nach) #-- z.B. iE_init()
        else
            #-- normale Kanten erzeugen
            Params = MakeParam(kk) 
            s = Symbol(typ,"_kante"); obj = getfield(FlexHyX, s)
            kante = obj(Param=Params, KL=knoten[von], KR=knoten[nach], Z=kk)    #-- z.B. iB_kante()
            push!(kanten, kante)
            append!(M, kante.M)
        end
    end

    M = sparse(diagm(M))

    #println("-------------------------------------------------------")

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
    end
    #-----------------------------------------------------

    #-- Erzeuge Inzidenzmatrix 
    IM, IP = inzidenz(knoten,kanten) 
    
    #-- Knoten bekommen Kanten zugewiesen
    for i in eachindex(knoten)
        in = findall(!iszero, IP[i,:])
        out = findall(!iszero, IM[i,:])
        if isempty(in) == false
            knoten[i].in = kanten[in]
        end
        if isempty(out) == false
            knoten[i].out = kanten[out]    
        end
    end

    #-- Erzeuge Zustandsvektor y und Indizes wo was steht in y
    y, idx_iflussL, idx_iflussR, idx_mflussL, idx_mflussR, idx_eflussL, idx_eflussR, idx_ele = netzwerk2y(knoten,kanten) 
    println("Anzahl Gleichungen: ",length(y))
    params = IM, IP, knoten, kanten, idx_iflussL, idx_iflussR, idx_mflussL, idx_mflussR, idx_eflussL, idx_eflussR, idx_ele

    #-- Jacobi Struktur
    dy = 0*y;
    if JPattern == 1
        jac_sparsity = Symbolics.jacobian_sparsity((dy, y) -> dgl!(dy, y, params, 0.0),copy(y), y)
        linsolver = KLUFactorization()
    else 
        jac_sparsity = nothing 
        linsolver = nothing
    end
    for i in eachindex(knoten)
        if hasfield(typeof(knoten[i].Param), :Jac_init) 
             knoten[i].Param.Jac_init = false
        end
    end
    for i in eachindex(kanten)
        if hasfield(typeof(kanten[i].Param), :Jac_init) 
            kanten[i].Param.Jac_init = false
        end
    end

    #-- konsistente AW berechnen 
    if AW == 1
        ind_alg = findall(x->x==0,M[diagind(M)]);
        dgl!(dy,y,params,0.0);
        println("Test Vorher: ",Base.maximum(abs.(dy[ind_alg])))
        y_alg = copy(y[ind_alg])
        g!(dy_alg,y_alg) = f_aw!(dy_alg,y_alg,ind_alg,y,params)
        res = nlsolve(g!,y_alg)
        y[ind_alg] = res.zero;
        dgl!(dy,y,params,0.0);
        println("Test Nachher: ",Base.maximum(abs.(dy[ind_alg])))
    end

    t0 = time()
    f = ODEFunction(dgl!; mass_matrix=M, jac_prototype = jac_sparsity)
    tspan = (0.0,simdauer)
    prob_ode = ODEProblem(f,y,tspan,params)

    if isempty(eventfile) 
        n_events = 0 
        sol = solve(prob_ode,Rodas5P(autodiff=true,diff_type=Val{:forward},linsolve=linsolver),progress=true, reltol=rtol,abstol=atol,dtmax=dtmax)
    else
        global n_events
        cb = VectorContinuousCallback(event_condition,event_affect!,n_events,affect_neg! = nothing)
        sol = solve(prob_ode,Rodas5P(autodiff=true,diff_type=Val{:forward},linsolve=linsolver), callback=cb, dense=false, progress=true, reltol=rtol, abstol=atol, dtmax=dtmax)
    end

    t1 = time()-t0
    println("CPU: ",t1)
    println(sol.retcode," nt=",size(sol.t))
    println(sol.destats)

    y = Leitsung_anhängen(sol,knoten,kanten,idx_iflussL,idx_iflussR,IM,IP)

    println("---------------- This was FlexHyX -----------------")
    return (y, sol.t, idx_ele)
end

function MakeParam(kk) 
    s = Symbol(kk["Typ"],"_Param"); P = getfield(FlexHyX, s)
    Param = P() #-- erzeuge Param z.B. iB_Param
    D = Dict()
    for ff in fieldnames(typeof(Param))
        if haskey(kk,String(ff))==true 
            D[String(ff)] = kk[String(ff)] #-- speichere geänderte Params in Dict
            #println(kk["Typ"]," --> ",String(ff),"=",kk[String(ff)])
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
        if (typeof(kanten[i]) <: Strom_Kante)
            idx_UL = kanten[i].KL.y.U
            idx_UR = kanten[i].KR.y.U
            if hasfield(typeof(kanten[i].y), :iL) == true
                idx_iL  = kanten[i].y.iL
                idx_iR  = kanten[i].y.iR
            else
                idx_iL  = kanten[i].y.i
                idx_iR  = idx_iL
            end
            UL = y[[idx_UL],:]; UR = y[[idx_UR],:]
            IL = y[[idx_iL],:]; IR = y[[idx_iR],:]
            P = UR.*IR - UL.*IL
            y = vcat(y,P)
        end
    end
    return y
end


