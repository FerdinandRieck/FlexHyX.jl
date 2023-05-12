function comp_jacstru(IP,IM,idx_ifluss, idx_mfluss, idx_efluss,elemente,neq)
    #Jacstru = sparse(zeros(Int,neq,neq))
    Jacstru = zeros(Int,neq,neq)
    k = 1; idx_k = Int[]; 
    for i = 1:length(elemente.knoten)
        idx_k = [idx_k; k]
        n_ele = length(elemente.knoten[i].M)
        n = k:k+n_ele-1
        if n[1] == n[end] n = k end
        if isdefined(elemente.knoten[i],:J)
            Jacstru[n,n] = Jacstru[n,n] + elemente.knoten[i].J
        end
        #-- Summe Flüsse 
        if isdefined(elemente.knoten[i],:J_fluss)
            J_fluss = elemente.knoten[i].J_fluss
            for ie = 1:length(elemente.knoten[i].J_fluss)
                if J_fluss[ie] == "sum_i"
                    iz = findall(!iszero, IP[i,:]); isR = idx_ifluss[iz,2]
                    iz = findall(!iszero, IM[i,:]); isL = idx_ifluss[iz,2]
                    Jacstru[k+ie-1,isR] .= 1; Jacstru[k+ie-1,isL] .= 1;
                end
                if J_fluss[ie] == "sum_m"
                    mz = findall(!iszero, IP[i,:]); msR = idx_mfluss[mz,2]
                    # !!! Fehler hier bei Netztwerk: Netz_Gas_Stromnetz.json !!!
                    mz = findall(!iszero, IM[i,:]); msL = idx_mfluss[mz,2] 
                    Jacstru[k+ie-1,msR] .= 1; Jacstru[k+ie-1,msL] .= 1;
                end
                if J_fluss[ie] == "sum_e"
                    ez = findall(!iszero, IP[i,:]); esR = idx_efluss[ez,2]
                    ez = findall(!iszero, IM[i,:]); esL = idx_efluss[ez,2]
                    Jacstru[k+ie-1,esR] .= 1; Jacstru[k+ie-1,esL] .= 1;
                end
            end
        end

        k = k+n_ele
    end

    for i = 1:length(elemente.kanten)
        #-- Linke und rechte Knoten auf eins setzen
        kk = elemente.kanten[i].Z
        KL = kk["VonNach"][1]; KR = kk["VonNach"][2]
        n_ele = length(elemente.kanten[i].M)

        n_KL = length(elemente.knoten[KL].M)
        J_KL = zeros(n_ele,n_KL); 
        ik = 0
        for ie = 1:n_ele
            for ff in fieldnames(typeof(elemente.knoten[KL].y))
                if ff != :Param
                    ik+=1
                    char_ff = first(String(ff))
                    idx_eq = String(Symbol("eq",ie))
                    vars = elemente.kanten[i].J_KL[idx_eq]
                    #if isnothing(findfirst(contains.(char_ff),vars))==false J_KL[ie,ik] = 1 end
                    if isempty(findall(x->x==char_ff,vars))==false J_KL[ie,ik] = 1 end
                    #!!! vielleicht auch mit haskey(collection, key) abfragen möglich?
                end
            end
            ik = 0
        end

        n_KR = length(elemente.knoten[KR].M)
        J_KR = zeros(n_ele,n_KR); 
        for ie = 1:n_ele
            for ff in fieldnames(typeof(elemente.knoten[KR].y))
                if ff != :Param
                    ik+=1
                    char_ff = first(String(ff))
                    idx_eq = String(Symbol("eq",ie))
                    vars = elemente.kanten[i].J_KR[idx_eq]
                    #if isnothing(findfirst(contains.(char_ff),vars))==false J_KR[ie,ik] = 1 end
                    if isempty(findall(x->x==char_ff,vars))==false J_KR[ie,ik] = 1 end
                end
            end
            ik = 0
        end

        n = k:k+n_ele-1
        n_KL = idx_k[KL]:idx_k[KL]+n_KL-1
        n_KR = idx_k[KR]:idx_k[KR]+n_KR-1
        if n[1] == n[end] n = n[1] end
        if n_KL[1] == n_KL[end] n_KL = n_KL[1] end
        if n_KR[1] == n_KR[end] n_KR = n_KR[1] end
        if length(J_KL)==1 
            Jacstru[n,n_KL] = Jacstru[n,n_KL] + J_KL[1]
        else
            Jacstru[n,n_KL] = Jacstru[n,n_KL] + J_KL
        end
        if length(J_KR)==1 
            Jacstru[n,n_KR] = Jacstru[n,n_KR] + J_KR[1]
        else
            Jacstru[n,n_KR] = Jacstru[n,n_KR] + J_KR
        end

        J_ele = elemente.kanten[i].J
        Jacstru[n,n] = Jacstru[n,n] + J_ele
        k = k+n_ele
    end
    return Jacstru
end

