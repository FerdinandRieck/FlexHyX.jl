global n_events = 4 # ???

function event_condition(out,y,t,integrator)
    IM, IP, elemente, i_fluss, m_fluss, e_fluss, idx_ifluss, idx_mfluss, idx_efluss, idx_ele, n_n, n_e = integrator.p
    #-- nur Nulldurchgänge von unten werden detektiert --!

    #-- Batterie --------------------------------------------------
    idx_ka5 = idx_ele["5i"][1]; idx_ka3 = idx_ele["3i"][1]
    io, status_laden = einaus_status(t,elemente.kanten[idx_ka5].Z["Schaltzeit"],0.1)
    io, status_entladen = einaus_status(t,elemente.kanten[idx_ka3].Z["Schaltzeit"],0.1)
    idx_ka1 = idx_ele["1qB"][1]
    idx_qB  = idx_ele["1qB"][2] 
    Batterie = elemente.kanten[idx_ka1];
    Batterie.Z["status_laden"] = status_laden
    Batterie.Z["status_entladen"] = status_entladen
    soc = y[idx_qB]/Batterie.Param.Q_max
    out[1] = 1
    #-- Event 1: SOC=0.4 wird unterschritten 
    if status_entladen == 1
        out[1] = 0.4-soc
    end
    #-- Event 2: Ueberladungsschutz Batterie 
    out[2] = soc-1.0
    #-- Event 3: SOC=0.8 wird unterschritten, Laden wieder einschalten 
    out[3] = 0.8-soc
    #-- Event 4: Batterien geladen, Entladen-Kante einschalten
    out[4] = soc-0.85

    #------ Elektrolyse --------------------------------------------------
    #idx_ka4 
end

function event_affect!(integrator, event_index)
    IM, IP, elemente, i_fluss, m_fluss, e_fluss, idx_ifluss, idx_mfluss, idx_efluss, idx_ele, n_n, n_e = integrator.p
    #--
    println("Event:",event_index,", t=",integrator.t)
    set_proposed_dt!(integrator,5.0) #-- neue Zeitschrittweite
    if event_index == 1
        idx_ka3 = idx_ele["3i"][1]; kk = elemente.kanten[idx_ka3].Z;
        kk["Schaltzeit"]=[-100.0,integrator.t+kk["Schaltdauer"]] #--!!!in Matlab t+1 #-- Batterie Entladen ausschalten
    end
    if event_index == 2
        idx_ka1 = idx_ele["1qB"][1]
        if elemente.kanten[idx_ka1].Z["status_laden"] == 1
            idx_ka5 = idx_ele["5i"][1]; kk = elemente.kanten[idx_ka5].Z;
            kk["Schaltzeit"]=[-100,integrator.t+kk["Schaltdauer"]]
        end
    end
    if event_index == 3
        idx_ka1 = idx_ele["1qB"][1]
        if elemente.kanten[idx_ka1].Z["status_laden"] == 0
            idx_ka5 = idx_ele["5i"][1]; kk = elemente.kanten[idx_ka5].Z;
            kk["Schaltzeit"]=[integrator.t+kk["Schaltdauer"]]
        end
    end
    if event_index == 4
        idx_ka1 = idx_ele["1qB"][1]
        if elemente.kanten[idx_ka1].Z["status_entladen"] == 0
            idx_ka3 = idx_ele["3i"][1]; kk = elemente.kanten[idx_ka3].Z;
            kk["Schaltzeit"]=[integrator.t+kk["Schaltdauer"]]
        end
    end

end
