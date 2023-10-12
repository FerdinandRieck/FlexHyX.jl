function dgl2!(dy,y,P,t) 
    IM, IP, knoten, kanten, idx_iflussL, idx_iflussR, idx_mflussL, idx_mflussR, idx_eflussL, idx_eflussR, idx_ele  = P
    
    array2netzwerk!_2(knoten,kanten,y)

    #-- jetzt alle Knoten und Kanten druchlaufen und Gleichungen erzeugen
    k = 1
    #-- Kantengleichungen
    for i in eachindex(kanten) 
        Kante!(dy,k,kanten[i],t) # Datentyp von kanten[i] bestimmt Kante! Funktion 
        n_ele = length(kanten[i].M)
        k = k+n_ele;
    end

    #-- Knotengleichungen
    for i in eachindex(knoten)
        if hasfield(typeof(knoten[i]), :sum_i)
            sum_i = IP[i,idx_iflussR[:,1]]'*y[idx_iflussR[:,2]] - IM[i,idx_iflussL[:,1]]'*y[idx_iflussL[:,2]] 
            knoten[i].sum_i = sum_i
        end
        if hasfield(typeof(knoten[i]), :sum_m)
            sum_m = IP[i,idx_mflussR[:,1]]'*y[idx_mflussR[:,2]] - IM[i,idx_mflussL[:,1]]'*y[idx_mflussL[:,2]]
            knoten[i].sum_m = sum_m
        end
        if hasfield(typeof(knoten[i]), :sum_e)
            sum_e = IP[i,idx_eflussR[:,1]]'*y[idx_eflussR[:,2]] - IM[i,idx_eflussL[:,1]]'*y[idx_eflussL[:,2]]
            knoten[i].sum_e = sum_e
        end
        Knoten!(dy,k,knoten[i],t) # Datentyp von knoten[i] bestimmt Knoten! Funktion
        n_ele = length(knoten[i].M)
        k = k+n_ele;
    end
end

function f_aw2!(dy_alg,y_alg,ind_alg,y,P)
    dy = 0*y
    y[ind_alg] = y_alg;
    dgl2!(dy,y,P,0.0)
    for i=1:length(y_alg) #-- keine Ahnung, warum das nicht mit dy_alg=dy[ind_alg] funktioniert
        dy_alg[i] = dy[ind_alg[i]];
    end
end
