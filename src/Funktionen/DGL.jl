function dgl!(dy,y,P,t) 
    IM, IP, elemente, sum_i, sum_m, sum_e, idx_iflussL, idx_iflussR, idx_mfluss, idx_efluss, idx_ele, n_n, n_e = P
    #-- Innerhalb von dgl-Funktion: alle benötigten Infos, Speicherplatz aus parameter vorhanden
    
    array2tuple!(elemente,y)    #??? Was passiert, wenn zeitschritt vom Solver verworfen wird??? Wird aktuelles y oder neues y verwendet?

    sum_i = IP[:,idx_iflussR[:,1]]*y[idx_iflussR[:,2]] - IM[:,idx_iflussL[:,1]]*y[idx_iflussL[:,2]] 
    sum_m = IP[:,idx_mfluss[:,1]]*y[idx_mfluss[:,2]] - IM[:,idx_mfluss[:,1]]*y[idx_mfluss[:,2]]
    sum_e = IP[:,idx_efluss[:,1]]*y[idx_efluss[:,2]] - IM[:,idx_efluss[:,1]]*y[idx_efluss[:,2]]

    #-- jetzt alle Knoten und Kanten druchlaufen und Gleichungen erzeugen

    k = 1
    #-- Knotengleichungen
    for i = 1:n_n 
        Knoten!(dy,k,sum_i[i],sum_m[i],sum_e[i],elemente.knoten[i]) # Datentyp von elemente.knoten[k] bestimmt Knoten! Funktion
        n_ele = length(elemente.knoten[i].M)
        k = k+n_ele;
    end

    #-- Kantengleichungen
    for i=1:n_e  
        Kante!(dy,k,elemente.kanten[i],t) # Datentyp von elemente.kanten[i] bestimmt Kante! Funktion 
        n_ele = length(elemente.kanten[i].M)
        k = k+n_ele;
    end
end

function f_aw!(dy_alg,y_alg,ind_alg,y,P)
    dy = 0*y
    y[ind_alg] = y_alg;
    dgl!(dy,y,P,0.0)
    for i=1:length(y_alg) #-- keine Ahnung, warum das nicht mit dy_alg=dy[ind_alg] funktioniert
        dy_alg[i] = dy[ind_alg[i]];
    end
end