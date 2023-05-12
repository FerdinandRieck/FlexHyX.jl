#=
global n_events = 1 # ???

function event_condition(out,y,t,integrator)
    IM, IP, elemente, i_fluss, m_fluss, e_fluss, idx_ifluss, idx_mfluss, idx_efluss, idx_ele, n_n, n_e = integrator.p
    #-- nur Nulldurchgänge von unten werden detektiert --!
    idx_ka = idx_ele["1qB"][1]
    idx_y  = idx_ele["1qB"][2] 
    Batterie = elemente.kanten[idx_ka];

    soc = y[idx_y]/Batterie.Param.Q_max
    out[1] = -soc+0.9
    #out[2] = 0.4-soc
end

function event_affect!(integrator, event_index)
    IM, IP, elemente, i_fluss, m_fluss, e_fluss, idx_ifluss, idx_mfluss, idx_efluss, idx_ele, n_n, n_e = integrator.p
    #--
    println("Event:",event_index,", t=",integrator.t)
    set_proposed_dt!(integrator,5.0) #-- neue Zeitschrittweite
    if event_index == 1
        idx_ka = idx_ele["5i"][1]; kk = elemente.kanten[idx_ka].Z;
        kk["Schaltzeit"]=[-100.0,integrator.t+kk["Schaltdauer"]] #-- Ladekante 5 ausschalten
    end
    #=
    if event_index == 2
        ii = findall(x->x==20,y2Nr)[1];  kk = kanten[y2ele[ii]];
        kk["Schaltzeit"]=[integrator.t+60.0] #-- Kante 20 BZ einschalten
        ii = findall(x->x==16,y2Nr)[1];  kk = kanten[y2ele[ii]];
        kk["Schaltzeit"]=[-100,integrator.t+60.0] #-- Kante 16 EL ausschalten
    end
    =#
end
=#