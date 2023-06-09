function dgl!(dy,y,P,t) 
    IM, IP, elemente, idx_iflussL, idx_iflussR, idx_mfluss, idx_efluss, idx_ele, n_n, n_e  = P
    
    array2netzwerk!(elemente,y)

    #-- jetzt alle Knoten und Kanten druchlaufen und Gleichungen erzeugen
    k = 1
    #-- Knotengleichungen
    for i = 1:n_n 
        if hasfield(typeof(elemente.knoten[i]), :sum_i)
            sum_i = IP[i,idx_iflussR[:,1]]'*y[idx_iflussR[:,2]] - IM[i,idx_iflussL[:,1]]'*y[idx_iflussL[:,2]] 
            elemente.knoten[i].sum_i = sum_i
        end
        if hasfield(typeof(elemente.knoten[i]), :sum_m)
            sum_m = IP[i,idx_mfluss[:,1]]'*y[idx_mfluss[:,2]] - IM[i,idx_mfluss[:,1]]'*y[idx_mfluss[:,2]]
            elemente.knoten[i].sum_m = sum_m
        end
        if hasfield(typeof(elemente.knoten[i]), :sum_e)
            sum_e = IP[i,idx_efluss[:,1]]'*y[idx_efluss[:,2]] - IM[i,idx_efluss[:,1]]'*y[idx_efluss[:,2]]
            elemente.knoten[i].sum_e = sum_e
        end
        Knoten!(dy,k,elemente.knoten[i],t) # Datentyp von elemente.knoten[i] bestimmt Knoten! Funktion
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
