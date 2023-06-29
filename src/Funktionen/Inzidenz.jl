function inzidenz(knoten,kanten)
    n_n = length(knoten);   n_e = length(kanten)
    IP = zeros(n_n,n_e); IM = zeros(n_n,n_e)
	for i = 1:n_e
        if haskey(kanten[i].Z,"RefKante")==true
            RefKante = kanten[i].Z["RefKante"]
            iv_RefKante = RefKante["VonNach"][1]
            in_RefKante = RefKante["VonNach"][2]
            IP[in_RefKante,i] = 1; IM[iv_RefKante,i] = 1;
        end
        iv = kanten[i].Z["VonNach"][1]; in = kanten[i].Z["VonNach"][2];
        IP[in,i] = 1; IM[iv,i] = 1;
    end
    return IM, IP
end