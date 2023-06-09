function inzidenz(knoten,kanten)
    n_n = size(knoten)[1];   n_e = size(kanten)[1]
    IP = zeros(n_n,n_e); IM = zeros(n_n,n_e)
	for i = 1:n_e
        if haskey(kanten[i],"RefKante")==true
            RefKante = kanten[i]["RefKante"]
            iv_RefKante = RefKante["VonNach"][1]
            in_RefKante = RefKante["VonNach"][2]
            IP[in_RefKante,i] = 1; IM[iv_RefKante,i] = 1;
        end
        iv = kanten[i]["VonNach"][1]; in = kanten[i]["VonNach"][2];
        IP[in,i] = 1; IM[iv,i] = 1;
    end
    return IM, IP
end