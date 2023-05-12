function inzidenz(knoten,kanten)
    n_n = size(knoten)[1];   n_e = size(kanten)[1]
    IP = zeros(n_n,n_e); IM = zeros(n_n,n_e)
	for i = 1:n_e
        iv = kanten[i]["VonNach"][1]; in = kanten[i]["VonNach"][2];;
        IP[in,i] = 1; IM[iv,i] = 1;
    end
    return IM, IP
end