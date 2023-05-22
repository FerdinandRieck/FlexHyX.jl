global n_events = 7

function event_condition(out,y,t,integrator)
    IM, IP, elemente, sum_i, sum_m, sum_e, idx_iflussL, idx_iflussR, idx_mfluss, idx_efluss, idx_ele, n_n, n_e = integrator.p
    
    #-- nur Nulldurchgänge von unten werden detektiert --!

    #-- Batterie --------------------------------------------------
    idx_ka1 = idx_ele["1q"][1]
    idx_qB  = idx_ele["1q"][2] 
    Batterie = elemente.kanten[idx_ka1]
    bat_soc = y[idx_qB]/Batterie.Param.Q_max
    batEmptyThres = 0.3
    batFullThres = 0.9

    #-- MHS --------------------------------------------------
    idx_kn15 = idx_ele["15Θ"][1]
    idx_Θ = idx_ele["15Θ"][2]
    MHS = elemente.knoten[idx_kn15]
    MHS_soc = y[idx_Θ]
    MHSEmptyThres = 0.1
    MHSFullThres = 0.9

    #-- Leistungsdifferenz -----------------------------------
    idx_U3 = idx_ele["3U"][2]
    idx_iPV = idx_ele["4i"][2]
    idx_iV = idx_ele["3i"][2]
    U3 = y[idx_U3]
    iPV = y[idx_iPV]
    P_gen = U3*iPV;
    iV = y[idx_iV];
    P_last = U3*iV;
    deltaP = P_gen-P_last;
    Batterie.Z["deltaPb"] = deltaP >= 0;
    #-----------------------------------

    ThresTol = 1e-3;

    if bat_soc <= batEmptyThres*(1 + ThresTol)
        Batterie.Z["batEmpty"] = true;
    else
        Batterie.Z["batEmpty"] = false;
    end
    if bat_soc >= batFullThres*(1 - ThresTol)
        Batterie.Z["batFull"] = true;
    else
        Batterie.Z["batFull"] = false;
    end
    
    
    if MHS_soc <= MHSEmptyThres*(1 + ThresTol)
        MHS.Z["MHSEmpty"] = true;
    else
        MHS.Z["MHSEmpty"] = false;
    end
    if MHS_soc >= MHSFullThres*(1 - ThresTol)
        MHS.Z["MHSFull"] = true;
    else
        MHS.Z["MHSFull"] = false;
    end
    

    #-----------------------------------
    #-- Event 1: Bat SOC=0.3 wird unterschritten
    out[1] = batEmptyThres - bat_soc

    #-- Event 2: Ueberladungsschutz Batterie, sollte nie erreicht werden
    out[2] = bat_soc - 1.0

    #-- Event 3: Batterien geladen
    out[3] = bat_soc - batFullThres

    #-- Event 4: deltaP wird positiv
    out[4] = deltaP

    #-- Event 5: deltaP wird negativ
    out[5] = -deltaP

    #-- Event 6: MHS leer
    out[6] = MHSEmptyThres - MHS_soc

    #-- Event 7: MHS voll
    out[7] = MHS_soc - MHSFullThres
    #-----------------------------------
end

function event_affect!(integrator, event_index)
    IM, IP, elemente, sum_i, sum_m, sum_e, idx_iflussL, idx_iflussR, idx_mfluss, idx_efluss, idx_ele, n_n, n_e  = integrator.p
    #--
    println("Event:",event_index,", t=",integrator.t)
    set_proposed_dt!(integrator,5.0) #-- neue Zeitschrittweite

    if event_index == 1
        regelung(idx_ele, elemente, integrator.t)
    end
    if event_index == 2
        println("State should not be reached!")
    end
    if event_index == 3
        regelung(idx_ele, elemente, integrator.t)
    end
    if event_index == 4
        idx_ka1 = idx_ele["1i"][1]
        Batterie = elemente.kanten[idx_ka1]
        Batterie.Z["deltaPb"] = true
        regelung(idx_ele, elemente, integrator.t)
    end
    if event_index == 5
        idx_ka1 = idx_ele["1i"][1]
        Batterie = elemente.kanten[idx_ka1]
        Batterie.Z["deltaPb"] = false
        regelung(idx_ele, elemente, integrator.t)
    end
    if event_index == 6
        regelung(idx_ele, elemente, integrator.t)
    end
    if event_index == 7
        regelung(idx_ele, elemente, integrator.t)
    end
end

function toggleSchalter(schalter, in, t)
    # An- und Ausschalten von Komponenten mit Pruefung, ob bereits an- oder aus
    x = false;
    if in == 1 || in == true
        x = true
    end
    
    state = einaus(t,schalter.Z["Schaltzeit"],0.1);
    on = state == 1
    off = state == 0
    
    if x && off # wenn angeschaltet werden soll und aktueller Zustand = aus
        schalter.Z["Schaltzeit"]=[t+0.1]
    end
    if ! x && on # wenn ausgeschaltet werden soll und aktueller Zustand = aus
        schalter.Z["Schaltzeit"]=[-100,t+0.1]
    end
end

function regelung(idx_ele, elemente, t)
    idx_ka1 = idx_ele["1i"][1]
    idx_ka7 = idx_ele["7i"][1]
    idx_ka9 = idx_ele["9i"][1]
    idx_ka18 = idx_ele["18i"][1]
    idx_kn15 = idx_ele["15Θ"][1]
    Batterie = elemente.kanten[idx_ka1]
    Grid_Schalter = elemente.kanten[idx_ka7]
    EL_Schalter = elemente.kanten[idx_ka9]
    BZ_Schalter = elemente.kanten[idx_ka18]
    MHS = elemente.knoten[idx_kn15]

    deltaPb = Batterie.Z["deltaPb"]
    batEmpty = Batterie.Z["batEmpty"]
    batFull = Batterie.Z["batFull"]
    MHSFull = MHS.Z["MHSFull"]
    MHSEmpty = MHS.Z["MHSEmpty"]

    if deltaPb
        if batFull
            if ! MHSFull
                toggleSchalter(EL_Schalter,1,t) #-- schalte AN
                toggleSchalter(BZ_Schalter,0,t) #-- schalte AUS
                toggleSchalter(Grid_Schalter,0,t) #-- schalte AUS
            else
                toggleSchalter(EL_Schalter,0,t) #-- schalte AUS
                toggleSchalter(BZ_Schalter,0,t) #-- schalte AUS
                toggleSchalter(Grid_Schalter,1,t) #-- schalte AN
            end
        else
            toggleSchalter(EL_Schalter,0,t) #-- schalte AUS
            toggleSchalter(BZ_Schalter,0,t) #-- schalte AUS
            toggleSchalter(Grid_Schalter,0,t) #-- schalte AUS
        end
    else
        if batEmpty
            if ! MHSEmpty
                toggleSchalter(EL_Schalter,0,t) #-- schalte AUS
                toggleSchalter(BZ_Schalter,1,t) #-- schalte AN
                toggleSchalter(Grid_Schalter,0,t) #-- schalte AUS
            else
                toggleSchalter(EL_Schalter,0,t) #-- schalte AUS
                toggleSchalter(BZ_Schalter,0,t) #-- schalte AUS
                toggleSchalter(Grid_Schalter,1,t) #-- schalte AN
            end
        else
            toggleSchalter(EL_Schalter,0,t) #-- schalte AUS
            toggleSchalter(BZ_Schalter,0,t) #-- schalte AUS
            toggleSchalter(Grid_Schalter,0,t) #-- schalte AUS
        end
    end

    println("Regelung aufgerufen");
end